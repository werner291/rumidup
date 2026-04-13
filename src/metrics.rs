use std::fmt;
use std::io::{Error as IoError, Write};

use serde::Serialize;

use crate::record::Flags;

/// Serializable metrics report. Produced from [`Metrics`] for JSON output.
/// Field names match the existing hand-rolled JSON and Picard-style TSV headers.
#[derive(Debug, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct MetricsReport<'a> {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<&'a str>,
    pub unpaired_reads_examined: usize,
    pub paired_reads_examined: usize,
    pub secondary_or_supplementary_rds: usize,
    pub unmapped_reads: usize,
    pub unpaired_read_duplicates: usize,
    pub unpaired_read_optical_duplicates: usize,
    pub read_pair_duplicates: usize,
    pub read_pair_optical_duplicates: usize,
    pub corrected_umis: usize,
    pub fraction_duplication: Option<f32>,
    pub estimated_library_size: u64,
}

/// Duplication metrics.
#[derive(Debug, Default)]
pub struct Metrics {
    unpaired_reads_examined: usize,
    paired_reads_examined: usize,
    secondary_or_supplementary_rds: usize,
    unmapped_reads: usize,
    unpaired_read_duplicates: usize,
    unpaired_read_optical_duplicates: usize,
    read_pair_duplicates: usize,
    read_pair_optical_duplicates: usize,
    corrected_umis: usize,
}

pub enum Status {
    UnpairedRead,
    PairedRead,
    SecondaryOrSupplementary,
    Unmapped,
    UnpairedDuplicate,
    UnpairedOpticalDuplicate,
    ReadpairDuplicate,
    ReadpairOpticalDuplicate,
    CorrectedUmi,
}

impl fmt::Display for Metrics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "UNPAIRED_READS_EXAMINED\tPAIRED_READS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tUNPAIRED_READ_OPTICAL_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tCORRECTED_UMIS\tFRACTION_DUPLICATION\tESTIMATED_LIBRARY_SIZE")?;
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.unpaired_reads_examined,
            self.paired_reads_examined,
            self.secondary_or_supplementary_rds,
            self.unmapped_reads,
            self.unpaired_read_duplicates,
            self.unpaired_read_optical_duplicates,
            self.read_pair_duplicates,
            self.read_pair_optical_duplicates,
            self.corrected_umis,
            self.fraction_duplication(),
            self.estimated_library_size()
        )
    }
}

#[inline]
fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

impl Metrics {
    pub fn fraction_duplication(&self) -> f32 {
        (self.read_pair_duplicates / 2 + self.unpaired_read_duplicates) as f32
            / (self.unpaired_reads_examined + self.paired_reads_examined / 2) as f32
    }

    /// Library size estimation ported from picard DuplicationMetrics
    /// https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/markduplicates/EstimateLibraryComplexity.java
    ///
    /// modified to count single end reads as well
    pub fn estimated_library_size(&self) -> u64 {
        let read_pairs = ((self.paired_reads_examined - self.read_pair_optical_duplicates) / 2
            + self.unpaired_reads_examined
            - self.unpaired_read_optical_duplicates) as f64;
        let unique_read_pairs = ((self.paired_reads_examined - self.read_pair_duplicates) / 2
            + self.unpaired_reads_examined
            - self.unpaired_read_duplicates) as f64;
        let read_pair_duplicates = read_pairs - unique_read_pairs;

        if read_pairs > 0.0 && read_pair_duplicates > 0.0 {
            let mut m = 1.0f64;
            let mut mm = 100.0f64;

            if unique_read_pairs >= read_pairs
                || f(m * unique_read_pairs, unique_read_pairs, read_pairs) < 0.0
            {
                panic!(
                    "Invalid values for pairs and unique pairs: {}, {}",
                    read_pairs, unique_read_pairs
                );
            }

            // find value of mm, large enough to act as other side for bisection method
            while f(mm * unique_read_pairs, unique_read_pairs, read_pairs) > 0.0 {
                mm *= 10.0;
            }

            // use bisection method (no more than 40 times) to find solution
            for _ in 0..40 {
                let r = (m + mm) / 2.0;
                let u = f(r * unique_read_pairs, unique_read_pairs, read_pairs);
                if u == 0.0 {
                    break;
                } else if u > 0.0 {
                    m = r;
                } else if u < 0.0 {
                    mm = r;
                }
            }

            (unique_read_pairs * (m + mm) / 2.0) as u64
        } else {
            0
        }
    }

    pub fn count_flags(&mut self, flags: Flags) {
        if flags.is_supplementary() || flags.is_secondary() {
            self.count(Status::SecondaryOrSupplementary);
        } else {
            if flags.is_unmapped() {
                self.count(Status::Unmapped);
            }

            if flags.is_segmented() {
                self.count(Status::PairedRead);
            } else {
                self.count(Status::UnpairedRead);
            }
        }
    }
    pub fn count_duplicate(&mut self, is_segmented: bool, is_optical: bool) {
        match (is_segmented, is_optical) {
            (true, true) => {
                self.count(Status::ReadpairOpticalDuplicate);
                self.count(Status::ReadpairDuplicate);
            }
            (true, false) => self.count(Status::ReadpairDuplicate),
            (false, true) => {
                self.count(Status::UnpairedDuplicate);
                self.count(Status::UnpairedOpticalDuplicate);
            }
            (false, false) => self.count(Status::UnpairedDuplicate),
        }
    }

    pub fn count(&mut self, status: Status) {
        self.count_many(status, 1);
    }

    pub fn count_many(&mut self, status: Status, count: usize) {
        match status {
            Status::UnpairedRead => self.unpaired_reads_examined += count,
            Status::PairedRead => self.paired_reads_examined += count,
            Status::SecondaryOrSupplementary => self.secondary_or_supplementary_rds += count,
            Status::Unmapped => self.unmapped_reads += count,
            Status::UnpairedDuplicate => self.unpaired_read_duplicates += count,
            Status::UnpairedOpticalDuplicate => self.unpaired_read_optical_duplicates += count,
            Status::ReadpairDuplicate => self.read_pair_duplicates += count,
            Status::ReadpairOpticalDuplicate => self.read_pair_optical_duplicates += count,
            Status::CorrectedUmi => self.corrected_umis += count,
        }
    }

    pub fn to_report<'a>(&self, sample_name: Option<&'a str>) -> MetricsReport<'a> {
        let frac = self.fraction_duplication();
        MetricsReport {
            sample_name,
            unpaired_reads_examined: self.unpaired_reads_examined,
            paired_reads_examined: self.paired_reads_examined,
            secondary_or_supplementary_rds: self.secondary_or_supplementary_rds,
            unmapped_reads: self.unmapped_reads,
            unpaired_read_duplicates: self.unpaired_read_duplicates,
            unpaired_read_optical_duplicates: self.unpaired_read_optical_duplicates,
            read_pair_duplicates: self.read_pair_duplicates,
            read_pair_optical_duplicates: self.read_pair_optical_duplicates,
            corrected_umis: self.corrected_umis,
            fraction_duplication: if frac.is_finite() { Some(frac) } else { None },
            estimated_library_size: self.estimated_library_size(),
        }
    }

    pub fn write_json<W: Write>(&self, w: W, sample_name: Option<&str>) -> Result<(), IoError> {
        serde_json::to_writer(w, &self.to_report(sample_name))
            .map_err(|e| IoError::new(std::io::ErrorKind::Other, e))
    }

    /// Legacy hand-rolled JSON, kept for testing equivalence with serde output.
    #[cfg(test)]
    fn write_json_legacy<W: Write>(&self, mut w: W) -> Result<(), IoError> {
        writeln!(w, "{{\"UNPAIRED_READS_EXAMINED\":{}, \"PAIRED_READS_EXAMINED\":{}, \"SECONDARY_OR_SUPPLEMENTARY_RDS\":{}, \"UNMAPPED_READS\":{}, \"UNPAIRED_READ_DUPLICATES\":{}, \"UNPAIRED_READ_OPTICAL_DUPLICATES\":{}, \"READ_PAIR_DUPLICATES\":{}, \"READ_PAIR_OPTICAL_DUPLICATES\":{}, \"CORRECTED_UMIS\":{}, \"FRACTION_DUPLICATION\":{}, \"ESTIMATED_LIBRARY_SIZE\":{}}}",
            self.unpaired_reads_examined,
            self.paired_reads_examined,
            self.secondary_or_supplementary_rds,
            self.unmapped_reads,
            self.unpaired_read_duplicates,
            self.unpaired_read_optical_duplicates,
            self.read_pair_duplicates,
            self.read_pair_optical_duplicates,
            self.corrected_umis,
            self.fraction_duplication(),
            self.estimated_library_size())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_metrics() -> Metrics {
        let mut m = Metrics::default();
        m.unpaired_reads_examined = 1000;
        m.paired_reads_examined = 5000;
        m.secondary_or_supplementary_rds = 50;
        m.unmapped_reads = 30;
        m.unpaired_read_duplicates = 100;
        m.unpaired_read_optical_duplicates = 10;
        m.read_pair_duplicates = 400;
        m.read_pair_optical_duplicates = 40;
        m.corrected_umis = 25;
        m
    }

    #[test]
    fn serde_matches_legacy_json_values() {
        let m = make_metrics();

        // Legacy output (may contain NaN for zero-read case, so test with real data)
        let mut legacy_buf = Vec::new();
        m.write_json_legacy(&mut legacy_buf).unwrap();
        let legacy: serde_json::Value = serde_json::from_slice(&legacy_buf).unwrap();

        // New serde output (without sample_name, to match legacy fields)
        let mut serde_buf = Vec::new();
        m.write_json(&mut serde_buf, None).unwrap();
        let new: serde_json::Value = serde_json::from_slice(&serde_buf).unwrap();

        // Compare every field that exists in the legacy output
        let legacy_obj = legacy.as_object().unwrap();
        let new_obj = new.as_object().unwrap();
        assert_eq!(
            legacy_obj.len(),
            new_obj.len(),
            "field count mismatch: legacy has {}, new has {}",
            legacy_obj.len(),
            new_obj.len()
        );
        for (key, legacy_val) in legacy_obj {
            let new_val = new_obj
                .get(key)
                .unwrap_or_else(|| panic!("missing key: {key}"));
            assert_eq!(
                legacy_val, new_val,
                "mismatch for {key}: legacy={legacy_val}, new={new_val}"
            );
        }
    }

    #[test]
    fn zero_reads_emits_null_fraction() {
        let m = Metrics::default();
        let mut buf = Vec::new();
        m.write_json(&mut buf, None).unwrap();
        let val: serde_json::Value = serde_json::from_slice(&buf).unwrap();
        assert!(val["FRACTION_DUPLICATION"].is_null());
    }

    #[test]
    fn json_includes_sample_name_when_provided() {
        let m = make_metrics();
        let mut buf = Vec::new();
        m.write_json(&mut buf, Some("sample_42")).unwrap();
        let val: serde_json::Value = serde_json::from_slice(&buf).unwrap();
        assert_eq!(val["SAMPLE_NAME"], "sample_42");
    }

    #[test]
    fn json_omits_sample_name_when_none() {
        let m = make_metrics();
        let mut buf = Vec::new();
        m.write_json(&mut buf, None).unwrap();
        let val: serde_json::Value = serde_json::from_slice(&buf).unwrap();
        assert!(val.get("SAMPLE_NAME").is_none());
    }
}

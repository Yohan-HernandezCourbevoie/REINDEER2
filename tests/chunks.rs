mod common;

use std::num::NonZero;

use common::{check_number_of_partitions, is_all_same, AutoRemoveDirectory, AutoRemoveFile};

use reindeer2::reindeer2::{read_fof_file, OutputFormat, Parameters, Reindeer2};

fn get_input_fof() -> String {
    String::from("tests/chunk_data/fof.txt")
}

fn get_output_dir() -> AutoRemoveDirectory {
    AutoRemoveDirectory::create_random()
}

fn get_query() -> String {
    String::from("tests/chunk_data/file1.fa")
}

#[test]
/// Checks that the number of chunk has no effect on RD2's output.
fn no_effect_chunk() {
    let input = get_input_fof();

    let (file_paths, nb_color) = read_fof_file(&input).unwrap();
    let mut parameters = Parameters {
        k: 31,
        m: 15,
        bf_size: 1u64 << 16,
        partition_number: 512,
        nb_color,
        abundance_number: NonZero::new(255).unwrap(),
        abundance_min: 0,
        abundance_max: NonZero::new(65024).unwrap(),
        dense_option: false,
        canonical: false,
        sampling_strategy: None,
    };
    assert_eq!(parameters.nb_color, 4);

    let partitions = check_number_of_partitions(
        &input,
        parameters.partition_number,
        parameters.abundance_number.get(),
        parameters.bf_size,
    );
    parameters.partition_number = partitions;

    let mut query_results = vec![];
    for chunks_size in 1..parameters.nb_color {
        let index_dir = get_output_dir();
        let mut index = Reindeer2::new(parameters.clone(), String::from(index_dir.filename()));
        let index_dir_name = index_dir.filename();
        let threshold = 0;
        let paths = file_paths.clone();

        index.build(paths, chunks_size, threshold).unwrap();

        let coverage = 0.5;
        let output_format = OutputFormat::AbundanceMatrix {
            format: reindeer2::reindeer2::MatrixFormat::Raw(None),
        };
        let query_output =
            AutoRemoveFile::create_from_path(format!("{}_query_results.csv", index_dir_name));
        let query_output_name = query_output.filename();

        index
            .query(
                &get_query(),
                index_dir_name,
                query_output_name,
                output_format,
                coverage,
            )
            .expect("Failed to query sequences");

        let data = std::fs::read_to_string(query_output_name).unwrap();
        query_results.push(data);
    }

    assert!(is_all_same(&query_results));
}

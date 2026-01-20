mod common;

use common::{check_number_of_partitions, is_all_same, AutoRemoveDirectory, AutoRemoveFile};

use reindeer2::reindeer2::{read_fof_file, OutputFormat, Reindeer2};

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
    let bf_size = 1u64 << 16;
    let partitions = 512;

    let kmer = 31;
    let minimizer = 15;
    let input = get_input_fof();
    let (file_paths, color_nb) = read_fof_file(&input).unwrap();
    let abundance = 255;
    let abundance_min = 0;
    let abundance_max = 65024;
    let dense_option = false;
    let canonical = false;

    let partitions = check_number_of_partitions(&input, partitions, abundance, bf_size);

    let index = Reindeer2::new(
        bf_size,
        partitions,
        kmer,
        minimizer,
        color_nb,
        abundance,
        abundance_min,
        abundance_max,
        dense_option,
        canonical,
    );

    let mut query_results = vec![];
    assert_eq!(color_nb, 4);
    for chunks_size in 1..color_nb {
        let index_dir = get_output_dir();
        let index_dir_name = index_dir.filename();
        let threshold = 0;
        let paths = file_paths.clone();

        index
            .build(paths, index_dir_name, chunks_size, threshold)
            .unwrap();

        let coverage = 0.5;
        let output_format = OutputFormat::AbundanceMatrix {
            format: reindeer2::reindeer2::AbundanceMatrixFormat::Raw(None),
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

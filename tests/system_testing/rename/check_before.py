import argparse


# INPUT = """index directory: "integration_test_index" parameters: Parameters { bf_size: 4294967296, partition_number: 510, k: 31, m: 15, nb_color: 2, abundance_number: 255, abundance_min: 0, abundance_max: 65024, dense_option: false, canonical: true, sampling_strategy: None, findere_z: 4, capacity: 2, } indexed filenames and k-mers: [ ("datasetB", 5850), ("datasetA", 5940), ]"""
EXPECTED_INFOS = [
    "abundance_max: 65024",
    "abundance_min: 0",
    "abundance_number: 255",
    "canonical: true",
    "capacity: 2",
    "dense_option: false",
    "findere_z: 4",
    "k: 31",
    "m: 15",
    "nb_color: 2",
    "partition_number: 510",
    "sampling_strategy: None",
    "bf_size: 4294967296",
]


def check(rd2_output: str):
    prefix = """index directory: "integration_test_index"\nparameters: Parameters {"""
    suffix = """,\n}\nindexed filenames and k-mers: [\n  ("datasetB", 5850),\n  ("datasetA", 5940),\n]\n\n"""

    assert rd2_output.startswith(prefix)
    rd2_output = rd2_output.removeprefix(prefix)
    assert rd2_output.endswith(suffix)
    rd2_output = rd2_output.removesuffix(suffix)

    infos = rd2_output.split(",")
    infos = [info.removeprefix("\n    ") for info in infos]
    infos.sort()
    expected = sorted(EXPECTED_INFOS)

    assert infos == expected


def read_content(filename: str):
    with open(filename, "r") as fichier:
        return fichier.read()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()

    filename = args.filename

    content = read_content(filename)
    check(content)


if __name__ == "__main__":
    main()

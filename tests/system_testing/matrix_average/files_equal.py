import argparse


def are_equal(first_file: str, second_file: str) -> bool:
    equal = True
    with open(first_file, "r") as f1, open(second_file, "r") as f2:
        for line1, line2 in zip(f1, f2):
            if line1 != line2:
                equal = False
                print(f"mismatch:\n{line2.strip()}\n{line1.strip()}")
    return equal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("first_file")
    parser.add_argument("second_file")
    args = parser.parse_args()

    first_file = args.first_file
    second_file = args.second_file

    return 0 if are_equal(first_file, second_file) else 1


if __name__ == "__main__":
    exit(main())

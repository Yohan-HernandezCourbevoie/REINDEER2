import argparse


def are_files_similar(file1: str, file2: str) -> bool:
    similar = True
    with open(file1, "r") as f1, open(file2, "r") as f2:
        first = True
        for line1, line2 in zip(f1, f2):
            if first:
                first = False
                if line1 != line2:
                    print(f"mismatch on first line:\n{line2.strip()}\n{line1.strip()}")
                    similar = False
                    break
            else:
                # first word must be the same, others words should constsis of same numbers
                # e.g. 1,5,8\t8,9 is similar to 1,8,5\t8,9
                words1, words2 = line1.split("\t"), line2.split("\t")
                if len(words1) != len(words2):
                    print(f"lenght mismatch:\n{words2}\n{words1}")
                    similar = False
                    break
                if words1[0] != words2[0]:
                    print(f"first word mismatch:\n{words2[0]}\n{words1[0]}")
                    similar = False
                    break
                words1 = words1[1:]
                words2 = words2[1:]
                for word1, word2 in zip(words1, words2):
                    word1 = sorted([int(x) for x in word1.split(",")])
                    word2 = sorted([int(x) for x in word2.split(",")])
                    if word1 != word2:
                        print(f"mismatch:\n{word2}\n{word1}")
                        similar = False
                        break
    return similar


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("first_file")
    parser.add_argument("second_file")
    args = parser.parse_args()

    first_file = args.first_file
    second_file = args.second_file

    return 0 if are_files_similar(first_file, second_file) else 1


if __name__ == "__main__":
    exit(main())

<!-- TODO this is a draft, needs polishing -->

# Public notions
These notions are manipulated by the users. 

## Chunk
A chunk is a group of files to be indexed. Files in a chunk are processed in parallel. The more files in a chunk, the faster the indexation (see "merge"), but the higher the RAM consumption.

## Merge
To allow for faster queries, REINDEER2 merges all the chunks into a single index at the end of the indexation step. The merge gets slower when the number of chunks increases.

# Private notions
User are not supposed to interact with these notions.

## Batch
A group of fasta records that are processed together.

## Saves
A save is a file written to disk after some indexation steps has been reached.
It allows for the continuation of the index after a crash, by simply starting at the last valid state reached.
Saves are written:
- after a chunk is done
- after a merge is done
Saves are removed after the indexation is done.
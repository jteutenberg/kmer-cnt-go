# Overview

This is a simple implementation of k-mer counting in golang similar to the Python, C++, and C implementations in by Heng Li (https://github.com/lh3/kmer-cnt).

This implementation is similar to kc-cpp1 in that counting is backed by the notoriously inefficient standard map, though it can scale with more cores as in kc-c4.

# Usage

Assuming you are using Heng Li's benchmark file: M_abscessus_HiSeq_10M.fa.gz

    mkdir $GOPATH/src/github.com/jteutenberg
    cd $GOPATH/src/github.com/jteutenberg
    git clone https://github.com/jteutenberg/kmer-cnt-go
    cd kmer-cnt-go
    go build kmer_count.go
    zcat M_abscessus_HiSeq_10M.fa.gz | ./kmer_count > kc-g1.out

# Performance

As noted above, this makes use of the standard golang map and this forms the bottleneck. By default it writes to 16 maps in parallel though on machines with more than 16 cores I expect this can usefully be increased to 32.

On an 8-core / 16 thread Zen machine this takes 33s, slightly less time than kc-c1.

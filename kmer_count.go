package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

//seqLineReader reads all lines from the input Reader, forwarding
//those that do not begin with '>' onto the lines channel
func seqLineReader(input io.Reader, lines chan<- string, done chan<- bool) {
	scanner := bufio.NewScanner(input)
	capacity := 4000*1000
	scanner.Buffer(make([]byte, capacity), capacity)
	for scanner.Scan() {
		line := scanner.Text()
		if line[0] != '>' {
			lines <- line
		}
	}
	done <- true
}

//splitNs reads strings from the lines channel, splits these by
//any non-ACGT characters, forwarding the tokens onto all the splitLines channels
func splitNs(minLength int, lines <-chan string, splitLines chan<- string, done chan<- bool) {
	for line, ok := <-lines; ok; line, ok = <-lines {
		start := 0
		for i,r := range line {
			if r != 'A' && r != 'C' && r != 'G' && r != 'T' {
				if start <= i-minLength {
					splitLines <- line[start:i]
				}
				start = i + 1
			}
		}
		//and send the last subsequence
		if start <= len(line)-minLength {
			splitLines <- line[start:]
		}
	}
	done <- true
}

//makeKmers converts strings of A,C,G,T into slices of uint64
//Each uint64 spans k characters and always represents the smaller of the k-mer and its
//reverse complement.
//This writes output k-mer lists to ALL channels provided in kmersOut
func makeKmers(k int, sequences <-chan string, kmersOut []chan []uint64, done chan<- bool) {
	toShift := uint((k-1)*2) //how far to shift a base to get to the reverse position
	mask := (^uint64(0)) >> uint(64 - k*2) //1s for all bits that will be used in a k-mer
	for seq, ok := <-sequences; ok; seq, ok = <-sequences {
		kmers := make([]uint64, len(seq)-k+1)
		var next uint64   //slide across bases as a k-mer
		var nextRC uint64 //the reverse-complement version
		for i, b := range []byte(seq) {
			//convert the character into a 0-3 value
			b = ((b >> 1) ^ ((b & 4) >> 2)) & 3
			bb := uint64(b)
			//shuffle it in to make the next k-mer and its reverse complement
			next = ((next << 2) | bb) & mask
			nextRC = ((nextRC >> 2) | ((^bb) << toShift) & mask) //XOR to complement
			if i < k-1 {
				//haven't seen enough bases to make the first k-mer yet
				continue
			}
			//select the minimum of the two k-mers
			minKmer := next
			if minKmer > nextRC {
				minKmer = nextRC
			}
			kmers[i-k+1] = minKmer
		}
		for _, kout := range kmersOut {
			kout <- kmers
		}
	}
	done <- true
}

//countKmers maintains a count of all k-mers in the input kmerSequences that matches filter when
//masked by filterMask
func countKmers(filter, filterMask uint64, counts map[uint64]int, kmerSequences <-chan []uint64, done chan<- bool) {
	for kmers, ok := <-kmerSequences; ok; kmers, ok = <-kmerSequences {
		for _, kmer := range kmers {
			if (kmer & filterMask) != filter {
				continue
			}
			//update the count
			if count, exists := counts[kmer]; exists {
				counts[kmer] = count + 1
			} else {
				counts[kmer] = 1
			}
		}
	}
	done <- true
}

//printHist writes a histogram of the number of entries of an int->int map with
//given value (to a maximum of 255) stdout.
func printHist(counts []map[uint64]int) {
	hist := make([]int, 256)
	for _, nextCounts := range counts {
		for _, count := range nextCounts {
			if count > 255 {
				count = 255
			}
			hist[count] += 1
		}
	}
	for i, count := range hist {
		if i > 0 {
			fmt.Println(i,"\t",count)
		}
	}
}

func main() {
	k := 31
	numMaps := 16 //only use powers of 2, or the map filtering will be wrong

	filterMask := uint64(numMaps-1) //split maps on using these lower bits

	// Pipeline stage 1: prepare an input stream from stdin
	doneLines := make(chan bool, 1)
	inputStrings := make(chan string, 3)
	go seqLineReader(os.Stdin, inputStrings, doneLines)

	// Pipeline stage 2: workers to split on and remove Ns etc.
	numWorkers := 4
	doneSplits := make(chan bool, numWorkers)
	sequences := make(chan string, 3)
	for i := 0; i < numWorkers; i++ {
		go splitNs(k, inputStrings, sequences, doneSplits)
	}

	// Pipeline stage 3: convert to kmers
	kmerWorkers := 8
	doneKmers := make(chan bool, kmerWorkers)
	kmerSequences := make([]chan []uint64, numMaps) //multiple output here
	for i := 0; i < numMaps; i++ {
		kmerSequences[i] = make(chan []uint64, 3)
	}
	for i := 0; i < kmerWorkers; i++ {
		go makeKmers(k, sequences, kmerSequences, doneKmers)
	}

	// Pipeline stage 4: count kmers from the good sequence strings
	allCounts := make([]map[uint64]int, numMaps)
	doneCounting := make(chan bool, numMaps)
	for i, seqs := range kmerSequences {
		allCounts[i] = make(map[uint64]int)
		go countKmers(uint64(i), filterMask, allCounts[i], seqs, doneCounting)
	}

	//let everything drain
	<-doneLines
	close(inputStrings)

	for i := 0; i < numWorkers; i++ {
		<-doneSplits
	}
	close(sequences)

	for i := 0; i < kmerWorkers; i++ {
		<-doneKmers
	}
	for _, s := range kmerSequences {
		close(s)
	}
	for i := 0; i < numMaps; i++ {
		<-doneCounting
	}

	//print the result
	printHist(allCounts)
}

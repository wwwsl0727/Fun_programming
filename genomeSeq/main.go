package main


import (
	"os"
	"log"
	"bufio"
	"fmt"
	"time"
)

func main() {
	/*
	m := make(map[string]int)
	m["ACT"] = 3
	m["GTGA"] = 6
  m["TA"] = 2
	m["adkfdjk "] = -4
	fmt.Println(MaxValue(m))
	*/
	/*
	text := "mamaliga"
	k := 2
	fmt.Println(FrequencyMap(text, k))
	*/
	/*
	text := "ACGTTGCATGTCGCATGATGCATGAGAGCT"
	k:=4
	fmt.Println(FrequentWords(text,k))
	*/
	/*
	genome := "ACGTACGT"
	k:=1
	L:=5
	t :=2
	fmt.Println(ClumpFinding(genome,k,L,t))
	*/
	/*
	genome := "CATGGGCATCGGCCATACGCC"
	fmt.Println(SkewArray(genome))
	*/
	/*
	genome := "AAAACGTCGAAAAA"
	fmt.Println(FrequencyMap(genome, 2))
	*/

	//Write all output to file 
	f, err := os.Create("/tmp/yourfile")
    check(err)
    defer f.Close()

    w := bufio.NewWriter(f)
	genomeFile:=os.Args[1]
	genome:=ReadGenomeFromFile(genomeFile)

	fmt.Println("We have read genome files")
	//fmt.Println("genome length",len(genome))
	k:=5
	L:=500
	t :=3

	//Apply the advanced clump finding algorithm
	outputAdvancedFile:=os.Args[2]
	patterns1:=AdvancedClumpFinding(genome,k,L,t)
	fmt.Println("number of patterns in advenced CF", len(patterns1))

	//Write results to file
	WritekmerToFile(patterns1, outputAdvancedFile)
	fmt.Println("We have written Advanced CF file")

	//Apply the original clump finding algorithm
	outputClumpFFile:=os.Args[3]
	patterns2:=ClumpFinding(genome,k,L,t)
	fmt.Println("number of patterns in original CF", len(patterns2))

	WritekmerToFile(patterns2, outputClumpFFile)
	fmt.Println("We have written original CF file")
	// ocfFileName:="ouput_CFkmers.txt"
	// acfFileName:="output_ACFkmers.txt"
	//
	// ocfkMers:=ReadkMersFromFile(ocfFileName)
	// acfkMers:=ReadkMersFromFile(acfFileName)
	// fmt.Println(SamePatterns(ocfkMers,acfkMers))
}

func SamePatterns(a,b []string) bool{
	//For each pattern in a check if there is a pair in b
	for _, pattern1 := range a{
		flag:=0
		for _,pattern2:= range b{
			//If there is one pattern in b matches a, check next pattern in a
			if pattern1==pattern2{
					flag=1
					break
				}
		}
		//If there is no pattern in b matches a, return false
		if flag ==0{
			return false
		}
	}
	return true
}

func ReadkMersFromFile(filename string) []string {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	kMerStrings := make([]string, 0)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		currentLine := scanner.Text()
		kMerStrings = append(kMerStrings, currentLine)
	}

	return kMerStrings
}

func ReadGenomeFromFile(filename string) string {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()


	fileinfo, err := file.Stat()
	if err != nil {
  	log.Fatal(err)
}

	var genome string

	fileSize := fileinfo.Size()
	reader := bufio.NewReaderSize(file, int(fileSize))
	genome,err = reader.ReadString('\n')

	return genome
}

func WritekmerToFile(patterns []string, filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Fatal("Cannot create file", err)
	}
	defer file.Close()

	for _, pattern := range patterns {
		fmt.Fprintln(file, pattern)
	}
}

func MinimumSkew(genome string) []int {
	minIndex := make([]int,0)
	skewArray := SkewArray(genome)
	minValue := MinValue(skewArray)
	for i:=range skewArray{
		if skewArray[i] ==minValue{
			minIndex = append(minIndex, i)
		}
	}
	return minIndex
}

func MinValue(array[]int)int{
	min := array[0]
	for i:=0;i<len(array);i++{
		if array[i] < min {
			min = array[i]
		}
	}
	return min
}
func SkewArray(genome string) []int {
	skewArray := make([]int,0)
	skewArray = append(skewArray, 0)
		for i := range(genome){
			switch genome[i]{
			case 'C':
						skewArray = append(skewArray, skewArray[i]-1)
			case 'G':
						skewArray = append(skewArray, skewArray[i]+1)
			case 'A':
						skewArray = append(skewArray, skewArray[i])
			 case 'T':
						skewArray = append(skewArray, skewArray[i])
			}
		}
		return skewArray
		}

//AdvancedClumpFinding: The only differenece between adjacent window is that the first-kmer of the former is deleted in the latter and the latter
//obtain a new window by adding a new character. Therefore, we don't need to calculate the FrequencyMap from the beginning.
//Instead, after sliding, we decrease freqMap[first-kmer] by 1 and add freqMap[lastKmer] by 1.  If the freqMap[lastKmer] >=0, add lastKmer to pattern.
func AdvancedClumpFinding(genome string, k, L, t int) []string {
	patterns := make([]string,0)
	n:=len(genome)
	var freqMap map[string]int

	start1 := time.Now()
	for i:=0;i<=n-L;i++{
		window:=genome[i:i+L]
		if i==0{
			freqMap = FrequencyMap(window, k)
			//Record the first window and add patterns whose appearances are bigger than k to list.
			for key,val := range freqMap{
				if val >=t{
					patterns=append(patterns,key)
				}
			}
		}else{
				//after sliding window, freqMap[firstkmer]--, freqMap[lastKmer]++, if it is t, then add it to the pattern
				//fmt.Println(i-1)
				firstKmer := genome[i-1:i-1+k]
				freqMap[firstKmer]--

				lastKmer :=window[L-k:L]
				freqMap[lastKmer]++
				if freqMap[lastKmer]==t{
					patterns =append(patterns, lastKmer)
				}
				}
		}
		elapsed1 := time.Since(start1)
		log.Printf("Advanced clump finding rule takes %s\n", elapsed1)

		start2 := time.Now()
		patterns = RemoveDuplicates(patterns)
		elapsed2 := time.Since(start2)
		log.Printf("Reduce duplicates in Advanced CF %s\n", elapsed2)
		return patterns
	}



func ClumpFinding(genome string, k, L, t int) []string {
	patterns := make([]string,0)
	n:=len(genome)
	start1:=time.Now()
	for i:=0;i<=n-L;i++{
		window:=genome[i:i+L]
		freqMap := FrequencyMap(window, k)
		for key,val := range freqMap{
			if val >=t{
				patterns=append(patterns,key)
			}
		}
	}
	elapsed1 := time.Since(start1)
	log.Printf("Original clump finding rule takes %s\n", elapsed1)

	start2 := time.Now()
	patterns = RemoveDuplicates(patterns)
	elapsed2 := time.Since(start2)
	log.Printf("Reduce duplicates in original CF %s\n", elapsed2)
	return patterns
}

func RemoveDuplicates(patterns []string) []string{
	distinctPatterns := make([]string,0)
	distinctPatternsDic := make(map[string] int)
	for i := range(patterns){
		//If the key firstly appears in the search
		key := patterns[i]
		_, ok := distinctPatternsDic[key]
		if ok==false{
				distinctPatternsDic[key] = 1
				distinctPatterns = append(distinctPatterns, key)
		}
		}
	return distinctPatterns
}

func MaxValue(freqmap map[string]int) int {
	max := 0
	first_time := true
	for _, v := range freqmap {
		if v > max || first_time == true {
			max = v
			first_time = false
		}
	}
	return max
}

//func KeyExistenceCheck(map map[string] int, key string) bool{

//}

func FrequencyMap(text string, k int) map[string]int{
	n:=len(text)
	dict :=make(map[string]int)

	for i:=0;i<=n-k;i++{
		key :=text[i:i+k]
		_, ok := dict[key]
		//If key exists in the key, value ++
		if ok  {
			dict[key] ++
		}else{
		//If key not exist in the key preciously, initialize the value as 1
			dict[key] = 1
		}
	}
	return dict
}

func FrequentWords(text string, k int) []string {
	dict := FrequencyMap(text, k)
	maxValue := MaxValue(dict)
	list := make([]string,0)
	for key,val := range dict{
		if 	val==maxValue{
			list = append(list,key)
		}
	}

	return list
}

Note to future self: make sure you read the bitwise flag the correct way

Number of reads in original file: 18,186,410
cat sorted_file.sam | grep -v "^@" | wc -l

Number of uniq reads: 13,719,048
cat C1_SE_uniqAlign_deduped.sam | grep -v "^@" | wc -l

Number of header lines: 65
cat C1_SE_uniqAlign_deduped.sam | grep -E "^@" | wc -l

Number of duplicates removed: 18,186,410 - 13,719,048 = 4,467,362

Number of wrong UMIs: 0

Uniq reads per chromosome: cat C1_SE_uniqAlign_deduped.sam | grep -v "^@" | cut -f 3 | uniq -c 
697508 1
564903 10
1220389 11
359951 12
467659 13
387239 14
437465 15
360923 16
517566 17
290506 18
571665 19
2787018 2
547615 3
589839 4
562160 5
510818 6
1113183 7
576463 8
627488 9
202002 MT
317853 X
2247 Y
3 JH584299.1
656 GL456233.2
6 GL456211.1
4 GL456221.1
1 GL456354.1
5 GL456210.1
4 GL456212.1
294 JH584304.1
2 GL456379.1
3 GL456367.1
1 GL456239.1
1 GL456383.1
5450 MU069435.1
1 GL456389.1
21 GL456370.1
1 GL456390.1
1 GL456382.1
17 GL456396.1
3 GL456368.1
3 MU069434.1
111 JH584295.1



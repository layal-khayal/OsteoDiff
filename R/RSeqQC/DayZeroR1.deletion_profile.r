pdf("/home/layal/SVN/DEXseq2Uniprot/results/DayZeroR1.deletion_profile.pdf")
pos=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100)
value=c(0,0,290,451,3472,2883,1277,1233,1390,1480,1345,1402,1509,1505,1737,1703,1847,1751,1864,1745,1810,1783,1841,1750,1880,1658,1708,1800,1785,1803,1798,1825,1797,1742,1906,1875,1746,1857,1806,1798,1704,1732,1783,1847,1803,1791,1729,1696,1696,1677,1645,1676,1744,1753,1703,1675,1859,1686,1732,1727,1777,1836,1786,1782,1697,1645,1645,1687,1853,1690,1688,1814,1766,1721,1683,1647,1786,1704,1745,1793,1799,1767,1868,1786,1765,1706,1778,1805,1645,1737,1675,1662,1695,1636,1701,1901,2443,2952,585,200,0)
plot(pos,value,type='b', col='blue',xlab="Read position (5'->3')", ylab='Deletion count')
dev.off()

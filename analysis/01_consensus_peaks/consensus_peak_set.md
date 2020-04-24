Purpose: Create consensus peak sets
===================================

To be conservative in our peak calls and be sure that we're dealing with replicatable peaks in the following analyses, here will will create a consensus peak .bed file for each DNA binding protein by taking only those peaks which overlap in all replicate experiments.

The resulting peaks will consist of peaks which overlapped by at least on base pair in each replicate and will use the `GenomicRanges::reduce` function to merge the peaks by taking the outer boundaries of overlapping peaks. This strategy may widen some peaks, but will ensure that each peak in the resulting peak set has evidence in all experiments performed for that DNA binding protein.

``` r
# For details on this function, see intersect_functions.R
# In this process we also filter to only the peaks on canonical chromosomes.
consensus_peaks <- create_consensus_peaks(broadpeakfilepath = "/Shares/rinn_class/data/k562_chip/results/bwa/mergedLibrary/macs/broadPeak")
```

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000253.2, GL383522.1, KI270706.1, KI270745.1, KI270751.1, KI270757.1, KI270761.1, KI270853.1, KI270861.1, KV766198.1, KZ208915.1, chrM
    ##   - in 'y': GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000221.1, GL000224.1, GL000225.1, GL000254.2, KI270438.1, KI270538.1, KI270723.1, KI270724.1, KI270726.1, KI270743.1, KI270816.1, KI270857.1, KN196487.1, ML143366.1, ML143367.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000253.2, GL383522.1, KI270706.1, KI270745.1, KI270751.1, KI270757.1, KI270761.1, KI270853.1, KI270861.1, KV766198.1, KZ208915.1, chrM
    ##   - in 'y': GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000221.1, GL000224.1, GL000225.1, GL000254.2, KI270438.1, KI270538.1, KI270723.1, KI270724.1, KI270726.1, KI270743.1, KI270816.1, KI270857.1, KN196487.1, ML143366.1, ML143367.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000253.2, GL000256.2, GL383522.1, KI270706.1, KI270745.1, KI270751.1, KI270757.1, KI270761.1, KI270853.1, KI270861.1, KI270908.1, KV766198.1, chrM, GL000214.1, GL000224.1, GL000225.1, KI270538.1, KN196487.1
    ##   - in 'y': KI270897.1, KI270905.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000253.2, GL000256.2, GL383522.1, KI270706.1, KI270745.1, KI270751.1, KI270757.1, KI270761.1, KI270853.1, KI270861.1, KI270908.1, KV766198.1, chrM, GL000214.1, GL000224.1, GL000225.1, KI270538.1, KN196487.1
    ##   - in 'y': KI270897.1, KI270905.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000216.2, GL000218.1, GL000221.1, GL000224.1, GL000225.1, GL000254.2, KI270438.1, KI270538.1, KI270723.1, KI270724.1, KI270726.1, KI270743.1, KI270816.1, KN196487.1, ML143366.1, ML143367.1
    ##   - in 'y': GL000208.1, GL000251.2, GL339449.2, KI270711.1, KI270728.1, KI270784.1, KI270810.1, KI270819.1, KI270822.1, KI270849.1, KI270879.1, KI270880.1, KN538364.1, KQ090026.1, KV575244.1, KV766197.1, KZ559111.1, ML143345.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000216.2, GL000218.1, GL000221.1, GL000224.1, GL000225.1, GL000254.2, KI270438.1, KI270538.1, KI270723.1, KI270724.1, KI270726.1, KI270743.1, KI270816.1, KN196487.1, ML143366.1, ML143367.1
    ##   - in 'y': GL000208.1, GL000251.2, GL339449.2, KI270711.1, KI270728.1, KI270784.1, KI270810.1, KI270819.1, KI270822.1, KI270849.1, KI270879.1, KI270880.1, KN538364.1, KQ090026.1, KV575244.1, KV766197.1, KZ559111.1, ML143345.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270754.1, KI270765.1, KV880764.1, ML143352.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000220.1, GL000225.1, KI270330.1, KI270442.1, KI270465.1, KI270538.1, KI270591.1, KI270709.1, KI270718.1, KI270722.1, KI270729.1, KI270736.1, KI270738.1, KI270742.1, KI270744.1, KI270751.1, KI270757.1, KN196487.1, KQ983257.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270754.1, KI270765.1, KV880764.1, ML143352.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000220.1, GL000225.1, KI270330.1, KI270442.1, KI270465.1, KI270538.1, KI270591.1, KI270709.1, KI270718.1, KI270722.1, KI270729.1, KI270736.1, KI270738.1, KI270742.1, KI270744.1, KI270751.1, KI270757.1, KN196487.1, KQ983257.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL339449.2, GL383522.1, GL383574.1, GL949746.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270583.1, KI270743.1, KI270765.1, KI270821.1, KI270871.1, KI270899.1, KI270904.1, KN538370.1, KN538372.1, KQ458384.1, KZ208907.1, KZ559109.1, ML143350.1, ML143352.1, ML143353.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000225.1, GL383563.3, GL877875.1, GL949752.1, JH159146.1, KI270726.1, KI270734.1, KI270762.1, KI270816.1, KI270819.1, KI270856.1, KI270868.1, KI270896.1, KI270903.1, KN196484.1, KQ031389.1, KQ090027.1, KV766193.1, KZ208922.1, KZ559100.1, KZ559103.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL339449.2, GL383522.1, GL383574.1, GL949746.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270583.1, KI270743.1, KI270765.1, KI270821.1, KI270871.1, KI270899.1, KI270904.1, KN538370.1, KN538372.1, KQ458384.1, KZ208907.1, KZ559109.1, ML143350.1, ML143352.1, ML143353.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000225.1, GL383563.3, GL877875.1, GL949752.1, JH159146.1, KI270726.1, KI270734.1, KI270762.1, KI270816.1, KI270819.1, KI270856.1, KI270868.1, KI270896.1, KI270903.1, KN196484.1, KQ031389.1, KQ090027.1, KV766193.1, KZ208922.1, KZ559100.1, KZ559103.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270442.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, KI270746.1, KI270751.1, KI270765.1, KI270880.1, KQ031384.1, KV766198.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': KI270713.1, KI270723.1, KI270908.1, KN196484.1, KZ208913.1, KZ208915.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270442.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, KI270746.1, KI270751.1, KI270765.1, KI270880.1, KQ031384.1, KV766198.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': KI270713.1, KI270723.1, KI270908.1, KN196484.1, KZ208913.1, KZ208915.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270765.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL000219.1, GL000255.2, GL000256.2, KI270709.1, KI270714.1, KI270731.1, KI270744.1, KI270844.1, KI270853.1, KI270879.1, KN196484.1, KQ090026.1, KQ458384.1, KV766198.1, KZ208912.1, KZ208915.1, KZ559112.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270765.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL000219.1, GL000255.2, GL000256.2, KI270709.1, KI270714.1, KI270731.1, KI270744.1, KI270844.1, KI270853.1, KI270879.1, KN196484.1, KQ090026.1, KQ458384.1, KV766198.1, KZ208912.1, KZ208915.1, KZ559112.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270765.1, GL000214.1, KI270709.1, KI270844.1, KI270853.1, KQ090026.1, KZ208912.1, chrM
    ##   - in 'y': GL000251.2, GL000253.2, KI270706.1, KI270707.1, KI270711.1, KI270742.1, KI270745.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KN538364.1, KZ208913.1, ML143345.1, ML143353.1, ML143371.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270765.1, GL000214.1, KI270709.1, KI270844.1, KI270853.1, KQ090026.1, KZ208912.1, chrM
    ##   - in 'y': GL000251.2, GL000253.2, KI270706.1, KI270707.1, KI270711.1, KI270742.1, KI270745.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KN538364.1, KZ208913.1, ML143345.1, ML143353.1, ML143371.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, GL000009.2, GL000218.1, GL000219.1, KI270731.1, KI270744.1, KI270844.1, KI270879.1, KN196484.1, KQ458384.1, KZ208912.1, KZ208915.1, chrM, KI270706.1, KI270707.1, KI270711.1, KI270742.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KZ208913.1, ML143371.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270712.1, KI270718.1, KI270728.1, KI270754.1, KI270856.1, KQ031389.1, KV880764.1, ML143352.1, ML143355.1, ML143365.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, GL000009.2, GL000218.1, GL000219.1, KI270731.1, KI270744.1, KI270844.1, KI270879.1, KN196484.1, KQ458384.1, KZ208912.1, KZ208915.1, chrM, KI270706.1, KI270707.1, KI270711.1, KI270742.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KZ208913.1, ML143371.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270712.1, KI270718.1, KI270728.1, KI270754.1, KI270856.1, KQ031389.1, KV880764.1, ML143352.1, ML143355.1, ML143365.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270765.1, GL000214.1, GL000218.1, GL000219.1, KI270709.1, KI270731.1, KI270844.1, KN196484.1, KQ090026.1, KQ458384.1, KZ208912.1, KZ208915.1, chrM, KI270706.1, KI270707.1, KI270711.1, KI270745.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KN538364.1, KZ208913.1, ML143345.1, ML143353.1, ML143371.1, ML143375.1, ML143377.1, ML143380.1, KI270712.1, KI270718.1, KI270754.1, KI270856.1, KQ031389.1, KV880764.1, ML143352.1, ML143355.1, ML143365.1, ML143379.1
    ##   - in 'y': GL383522.1, KI270337.1, KI270467.1, KV575244.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270765.1, GL000214.1, GL000218.1, GL000219.1, KI270709.1, KI270731.1, KI270844.1, KN196484.1, KQ090026.1, KQ458384.1, KZ208912.1, KZ208915.1, chrM, KI270706.1, KI270707.1, KI270711.1, KI270745.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270905.1, KN538364.1, KZ208913.1, ML143345.1, ML143353.1, ML143371.1, ML143375.1, ML143377.1, ML143380.1, KI270712.1, KI270718.1, KI270754.1, KI270856.1, KQ031389.1, KV880764.1, ML143352.1, ML143355.1, ML143365.1, ML143379.1
    ##   - in 'y': GL383522.1, KI270337.1, KI270467.1, KV575244.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270713.1, KI270714.1, KI270742.1, KI270765.1, KI270832.1, KN538370.1, KQ031389.1, KV880764.1, ML143345.1, ML143355.1, ML143372.1, ML143379.1
    ##   - in 'y': GL000256.2, GL383520.2, KI270442.1, KI270466.1, KI270467.1, KI270728.1, KI270850.1, ML143366.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270713.1, KI270714.1, KI270742.1, KI270765.1, KI270832.1, KN538370.1, KQ031389.1, KV880764.1, ML143345.1, ML143355.1, ML143372.1, ML143379.1
    ##   - in 'y': GL000256.2, GL383520.2, KI270442.1, KI270466.1, KI270467.1, KI270728.1, KI270850.1, ML143366.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270783.1, KQ090017.1, KV575244.1
    ##   - in 'y': GL000205.2, GL000216.2, GL000225.1, GL000256.2, GL000258.2, GL383542.1, GL949746.1, JH159137.1, KI270706.1, KI270721.1, KI270733.1, KI270734.1, KI270743.1, KI270749.1, KI270750.1, KI270818.1, KI270821.1, KI270822.1, KI270827.1, KI270831.1, KI270849.1, KI270862.1, KI270869.1, KI270879.1, KI270881.1, KI270896.1, KI270904.1, KI270925.1, KN196478.1, KN196479.1, KN196482.1, KN538364.1, KQ090023.1, KV766198.1, KZ208912.1, KZ208913.1, KZ559103.1, KZ559109.1, ML143366.1, ML143367.1, ML143377.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270783.1, KQ090017.1, KV575244.1
    ##   - in 'y': GL000205.2, GL000216.2, GL000225.1, GL000256.2, GL000258.2, GL383542.1, GL949746.1, JH159137.1, KI270706.1, KI270721.1, KI270733.1, KI270734.1, KI270743.1, KI270749.1, KI270750.1, KI270818.1, KI270821.1, KI270822.1, KI270827.1, KI270831.1, KI270849.1, KI270862.1, KI270869.1, KI270879.1, KI270881.1, KI270896.1, KI270904.1, KI270925.1, KN196478.1, KN196479.1, KN196482.1, KN538364.1, KQ090023.1, KV766198.1, KZ208912.1, KZ208913.1, KZ559103.1, KZ559109.1, ML143366.1, ML143367.1, ML143377.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000252.2, GL339449.2, GL383520.2, GL383556.1, GL383578.2, KI270731.1, KI270783.1, KI270821.1, KI270849.1, KI270924.1, KN196478.1, KN538364.1, KQ031389.1, KQ090023.1, KV880764.1, KZ208907.1, KZ559105.1, ML143380.1, chrM
    ##   - in 'y': GL000008.2, GL000218.1, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270438.1, KI270442.1, KI270538.1, KI270723.1, KI270743.1, KI270751.1, KI270770.1, KI270832.1, KN196487.1, KN538372.1, ML143344.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000252.2, GL339449.2, GL383520.2, GL383556.1, GL383578.2, KI270731.1, KI270783.1, KI270821.1, KI270849.1, KI270924.1, KN196478.1, KN538364.1, KQ031389.1, KQ090023.1, KV880764.1, KZ208907.1, KZ559105.1, ML143380.1, chrM
    ##   - in 'y': GL000008.2, GL000218.1, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270438.1, KI270442.1, KI270538.1, KI270723.1, KI270743.1, KI270751.1, KI270770.1, KI270832.1, KN196487.1, KN538372.1, ML143344.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000251.2, GL000252.2, GL000258.2, GL383520.2, GL383563.3, GL383578.2, KI270711.1, KI270714.1, KI270721.1, KI270722.1, KI270724.1, KI270731.1, KI270741.1, KI270754.1, KI270782.1, KI270783.1, KI270787.1, KI270803.1, KI270819.1, KI270827.1, KI270830.1, KI270850.1, KI270894.1, KI270899.1, KI270905.1, KI270924.1, KN196478.1, KN196484.1, KN538364.1, KQ031389.1, KQ090023.1, KQ458384.1, KQ458385.1, KQ983257.1, KV880764.1, KV880768.1, KZ208907.1, KZ559105.1, ML143345.1, ML143373.1, GL000008.2, GL000218.1, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270538.1, KI270723.1, KI270743.1, KI270770.1, KI270832.1, KN196487.1, KN538372.1, ML143344.1, chrY
    ##   - in 'y': KI270333.1, KI270337.1, KI270466.1, KI270467.1, KI270737.1, KI270822.1, KI270896.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000251.2, GL000252.2, GL000258.2, GL383520.2, GL383563.3, GL383578.2, KI270711.1, KI270714.1, KI270721.1, KI270722.1, KI270724.1, KI270731.1, KI270741.1, KI270754.1, KI270782.1, KI270783.1, KI270787.1, KI270803.1, KI270819.1, KI270827.1, KI270830.1, KI270850.1, KI270894.1, KI270899.1, KI270905.1, KI270924.1, KN196478.1, KN196484.1, KN538364.1, KQ031389.1, KQ090023.1, KQ458384.1, KQ458385.1, KQ983257.1, KV880764.1, KV880768.1, KZ208907.1, KZ559105.1, ML143345.1, ML143373.1, GL000008.2, GL000218.1, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270538.1, KI270723.1, KI270743.1, KI270770.1, KI270832.1, KN196487.1, KN538372.1, ML143344.1, chrY
    ##   - in 'y': KI270333.1, KI270337.1, KI270466.1, KI270467.1, KI270737.1, KI270822.1, KI270896.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000214.1, GL000252.2, GL000253.2, GL000258.2, GL339449.2, GL383520.2, GL383556.1, GL383563.3, GL383578.2, KI270711.1, KI270714.1, KI270720.1, KI270721.1, KI270722.1, KI270724.1, KI270731.1, KI270741.1, KI270750.1, KI270754.1, KI270782.1, KI270783.1, KI270787.1, KI270802.1, KI270803.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270850.1, KI270853.1, KI270894.1, KI270899.1, KI270905.1, KI270924.1, KN196478.1, KN196484.1, KN538364.1, KQ031389.1, KQ090023.1, KQ458384.1, KQ458385.1, KV766198.1, KV880764.1, KZ208907.1, KZ559105.1, ML143345.1, ML143373.1, GL000008.2, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270442.1, KI270538.1, KI270723.1, KI270743.1, KI270770.1, KI270832.1, ML143344.1, chrY, KI270737.1, KI270822.1, KI270896.1
    ##   - in 'y': GL000194.1, KI270746.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000214.1, GL000252.2, GL000253.2, GL000258.2, GL339449.2, GL383520.2, GL383556.1, GL383563.3, GL383578.2, KI270711.1, KI270714.1, KI270720.1, KI270721.1, KI270722.1, KI270724.1, KI270731.1, KI270741.1, KI270750.1, KI270754.1, KI270782.1, KI270783.1, KI270787.1, KI270802.1, KI270803.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270850.1, KI270853.1, KI270894.1, KI270899.1, KI270905.1, KI270924.1, KN196478.1, KN196484.1, KN538364.1, KQ031389.1, KQ090023.1, KQ458384.1, KQ458385.1, KV766198.1, KV880764.1, KZ208907.1, KZ559105.1, ML143345.1, ML143373.1, GL000008.2, GL000220.1, GL000221.1, KI270310.1, KI270315.1, KI270442.1, KI270538.1, KI270723.1, KI270743.1, KI270770.1, KI270832.1, ML143344.1, chrY, KI270737.1, KI270822.1, KI270896.1
    ##   - in 'y': GL000194.1, KI270746.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, KN538364.1, KQ458384.1, KZ559105.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000257.2, GL383530.1, KI270707.1, KI270722.1, KI270723.1, KI270724.1, KI270750.1, KI270830.1, KI270878.1, KV880763.1, KV880768.1, ML143378.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, KN538364.1, KQ458384.1, KZ559105.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000257.2, GL383530.1, KI270707.1, KI270722.1, KI270723.1, KI270724.1, KI270750.1, KI270830.1, KI270878.1, KV880763.1, KV880768.1, ML143378.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, GL000225.1, GL000251.2, KI270315.1, KI270538.1, KI270720.1, KI270726.1, KI270728.1, KI270734.1, KI270761.1, KI270802.1, KI270803.1, KI270819.1, KI270831.1, KI270832.1, KI270907.1, KN538372.1, KQ983257.1, KV880764.1, KZ208912.1, ML143366.1, ML143373.1
    ##   - in 'y': GL000221.1, GL000226.1, GL383578.2, GL383580.2, KI270310.1, KI270707.1, KI270733.1, KI270787.1, KI270868.1, KI270880.1, KN196484.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, GL000225.1, GL000251.2, KI270315.1, KI270538.1, KI270720.1, KI270726.1, KI270728.1, KI270734.1, KI270761.1, KI270802.1, KI270803.1, KI270819.1, KI270831.1, KI270832.1, KI270907.1, KN538372.1, KQ983257.1, KV880764.1, KZ208912.1, ML143366.1, ML143373.1
    ##   - in 'y': GL000221.1, GL000226.1, GL383578.2, GL383580.2, KI270310.1, KI270707.1, KI270733.1, KI270787.1, KI270868.1, KI270880.1, KN196484.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL000251.2, KI270442.1, KI270712.1, KI270713.1, KI270732.1, KI270871.1, KV766198.1, ML143352.1, ML143355.1, ML143366.1, chr10, chr13, chr18, chr19, chr2, chr22, chr4, chrX
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000225.1, KI270590.1, KI270742.1, KI270744.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL000251.2, KI270442.1, KI270712.1, KI270713.1, KI270732.1, KI270871.1, KV766198.1, ML143352.1, ML143355.1, ML143366.1, chr10, chr13, chr18, chr19, chr2, chr22, chr4, chrX
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000225.1, KI270590.1, KI270742.1, KI270744.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, KI270712.1, KI270713.1, KI270732.1, KI270871.1, KV766198.1, ML143355.1, ML143366.1, GL000214.1, GL000225.1, KI270590.1, KI270742.1, KN196487.1
    ##   - in 'y': KI270754.1, KI270869.1, ML143353.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, KI270712.1, KI270713.1, KI270732.1, KI270871.1, KV766198.1, ML143355.1, ML143366.1, GL000214.1, GL000225.1, KI270590.1, KI270742.1, KN196487.1
    ##   - in 'y': KI270754.1, KI270869.1, ML143353.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL000251.2, KI270538.1, KI270712.1, KI270732.1, KI270757.1, KI270871.1, KV766198.1, ML143366.1, GL000225.1, KI270590.1, KN196487.1, KI270869.1
    ##   - in 'y': KI270718.1, KI270728.1, KI270738.1, KI270782.1, KI270899.1, KN538364.1, KN538370.1, KQ031389.1, KV880764.1, ML143345.1, ML143370.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL000251.2, KI270538.1, KI270712.1, KI270732.1, KI270757.1, KI270871.1, KV766198.1, ML143366.1, GL000225.1, KI270590.1, KN196487.1, KI270869.1
    ##   - in 'y': KI270718.1, KI270728.1, KI270738.1, KI270782.1, KI270899.1, KN538364.1, KN538370.1, KQ031389.1, KV880764.1, ML143345.1, ML143370.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383556.1, GL383563.3, GL383575.2, GL877875.1, GL949746.1, KI270582.1, KI270711.1, KI270718.1, KI270735.1, KI270787.1, KI270821.1, KI270836.1, KI270845.1, KI270871.1, KI270924.1, KI270934.1, KQ458384.1, KZ208913.1, ML143366.1
    ##   - in 'y': GL000216.2, GL000224.1, GL000254.2, GL383555.2, GL383581.2, JH636055.2, KI270438.1, KI270709.1, KI270721.1, KI270742.1, KI270754.1, KI270762.1, KI270772.1, KI270783.1, KI270792.1, KI270846.1, KI270854.1, KI270869.1, KI270879.1, KI270903.1, KI270904.1, KI270937.1, KN196479.1, KQ458383.1, KQ983257.1, KZ208922.1, KZ559111.1, ML143353.1, ML143358.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383556.1, GL383563.3, GL383575.2, GL877875.1, GL949746.1, KI270582.1, KI270711.1, KI270718.1, KI270735.1, KI270787.1, KI270821.1, KI270836.1, KI270845.1, KI270871.1, KI270924.1, KI270934.1, KQ458384.1, KZ208913.1, ML143366.1
    ##   - in 'y': GL000216.2, GL000224.1, GL000254.2, GL383555.2, GL383581.2, JH636055.2, KI270438.1, KI270709.1, KI270721.1, KI270742.1, KI270754.1, KI270762.1, KI270772.1, KI270783.1, KI270792.1, KI270846.1, KI270854.1, KI270869.1, KI270879.1, KI270903.1, KI270904.1, KI270937.1, KN196479.1, KQ458383.1, KQ983257.1, KZ208922.1, KZ559111.1, ML143353.1, ML143358.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270731.1, KI270737.1, KI270765.1, KI270841.1, ML143352.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000221.1, GL000226.1, KI270384.1, KI270465.1, KI270510.1, KI270515.1, KI270517.1, KI270519.1, KI270538.1, KI270590.1, KI270724.1, KI270736.1, KI270746.1, KI270747.1, KI270750.1, KI270751.1, KI270756.1, KI270772.1, KI270846.1, KN196483.1, KQ090024.1, ML143357.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270731.1, KI270737.1, KI270765.1, KI270841.1, ML143352.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000221.1, GL000226.1, KI270384.1, KI270465.1, KI270510.1, KI270515.1, KI270517.1, KI270519.1, KI270538.1, KI270590.1, KI270724.1, KI270736.1, KI270746.1, KI270747.1, KI270750.1, KI270751.1, KI270756.1, KI270772.1, KI270846.1, KN196483.1, KQ090024.1, ML143357.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270442.1, KI270709.1, KI270730.1, KI270734.1, KI270756.1, KI270783.1, KI270784.1, KI270857.1, KN538360.1, KN538364.1, KN538370.1, KV575244.1, KV880764.1, KV880768.1, ML143345.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270714.1, KI270782.1, KI270830.1, KI270849.1, KI270856.1, KI270861.1, KI270879.1, KI270880.1, KI270905.1, KN196484.1, KZ208912.1, ML143366.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270442.1, KI270709.1, KI270730.1, KI270734.1, KI270756.1, KI270783.1, KI270784.1, KI270857.1, KN538360.1, KN538364.1, KN538370.1, KV575244.1, KV880764.1, KV880768.1, ML143345.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270714.1, KI270782.1, KI270830.1, KI270849.1, KI270856.1, KI270861.1, KI270879.1, KI270880.1, KI270905.1, KN196484.1, KZ208912.1, ML143366.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000258.2, KI270466.1, KI270467.1, KI270830.1, KZ208908.1, KZ208912.1, KZ559109.1, ML143352.1, ML143353.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000252.2, GL383533.1, GL383563.3, JH159136.1, KI270330.1, KI270442.1, KI270710.1, KI270726.1, KI270731.1, KI270732.1, KI270733.1, KI270735.1, KI270738.1, KI270745.1, KI270753.1, KI270763.1, KI270782.1, KI270799.1, KI270803.1, KI270816.1, KI270856.1, KI270875.1, KI270878.1, KI270896.1, KI270900.1, KI270904.1, KI270925.1, KN196484.1, KQ031389.1, KQ090022.1, KQ090026.1, KQ759759.1, KV766193.1, KV880765.1, KZ208922.1, ML143364.1, ML143367.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000258.2, KI270466.1, KI270467.1, KI270830.1, KZ208908.1, KZ208912.1, KZ559109.1, ML143352.1, ML143353.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, GL000252.2, GL383533.1, GL383563.3, JH159136.1, KI270330.1, KI270442.1, KI270710.1, KI270726.1, KI270731.1, KI270732.1, KI270733.1, KI270735.1, KI270738.1, KI270745.1, KI270753.1, KI270763.1, KI270782.1, KI270799.1, KI270803.1, KI270816.1, KI270856.1, KI270875.1, KI270878.1, KI270896.1, KI270900.1, KI270904.1, KI270925.1, KN196484.1, KQ031389.1, KQ090022.1, KQ090026.1, KQ759759.1, KV766193.1, KV880765.1, KZ208922.1, ML143364.1, ML143367.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, KI270330.1, KI270707.1, KI270712.1, KI270714.1, KI270784.1, KI270849.1, KI270876.1, KQ458385.1, KZ208912.1, KZ208913.1, ML143344.1, ML143352.1
    ##   - in 'y': GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000254.2, KI270709.1, KI270711.1, KI270733.1, KI270743.1, KI270751.1, KI270821.1, KI270827.1, KI270832.1, KI270872.1, KQ983257.1, KV575244.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, KI270330.1, KI270707.1, KI270712.1, KI270714.1, KI270784.1, KI270849.1, KI270876.1, KQ458385.1, KZ208912.1, KZ208913.1, ML143344.1, ML143352.1
    ##   - in 'y': GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000252.2, GL000254.2, KI270709.1, KI270711.1, KI270733.1, KI270743.1, KI270751.1, KI270821.1, KI270827.1, KI270832.1, KI270872.1, KQ983257.1, KV575244.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000258.2, GL383578.2, GL877875.1, KI270712.1, KI270733.1, KI270765.1, KI270782.1, KI270810.1, KI270818.1, KI270880.1, KN538364.1, KZ559111.1, ML143352.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000220.1, GL000224.1, KI270709.1, KI270831.1, KN196484.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000258.2, GL383578.2, GL877875.1, KI270712.1, KI270733.1, KI270765.1, KI270782.1, KI270810.1, KI270818.1, KI270880.1, KN538364.1, KZ559111.1, ML143352.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000220.1, GL000224.1, KI270709.1, KI270831.1, KN196484.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000253.2, GL383526.1, GL383567.1, GL877875.1, GL949746.1, KI270706.1, KI270721.1, KI270728.1, KI270783.1, KI270784.1, KI270813.1, KI270822.1, KI270844.1, KI270853.1, KI270879.1, KI270899.1, KI270903.1, KI270904.1, KI270908.1, KN196478.1, KQ090026.1, KV766198.1, KV880763.1, KZ559103.1, KZ559112.1, ML143345.1, ML143352.1
    ##   - in 'y': GL000216.2, GL000218.1, GL000219.1, GL000220.1, KI270438.1, KI270723.1, KI270733.1, KI270754.1, KN196484.1, KN196487.1, KV575244.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000253.2, GL383526.1, GL383567.1, GL877875.1, GL949746.1, KI270706.1, KI270721.1, KI270728.1, KI270783.1, KI270784.1, KI270813.1, KI270822.1, KI270844.1, KI270853.1, KI270879.1, KI270899.1, KI270903.1, KI270904.1, KI270908.1, KN196478.1, KQ090026.1, KV766198.1, KV880763.1, KZ559103.1, KZ559112.1, ML143345.1, ML143352.1
    ##   - in 'y': GL000216.2, GL000218.1, GL000219.1, GL000220.1, KI270438.1, KI270723.1, KI270733.1, KI270754.1, KN196484.1, KN196487.1, KV575244.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000251.2, GL000256.2, GL339449.2, KI270467.1, KI270832.1, KI270857.1, KI270899.1, KV575245.1, ML143352.1, ML143361.1
    ##   - in 'y': KI270709.1, KI270853.1, KQ983257.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000251.2, GL000256.2, GL339449.2, KI270467.1, KI270832.1, KI270857.1, KI270899.1, KV575245.1, ML143352.1, ML143361.1
    ##   - in 'y': KI270709.1, KI270853.1, KQ983257.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270744.1, KI270745.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270821.1, KI270871.1, KI270905.1, KZ208915.1, KZ559112.1, ML143371.1, chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000220.1, GL000224.1, KI270438.1, KI270709.1, KI270719.1, KI270733.1, KI270738.1, KI270754.1, KI270757.1, KI270908.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270744.1, KI270745.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270821.1, KI270871.1, KI270905.1, KZ208915.1, KZ559112.1, ML143371.1, chrM
    ##   - in 'y': GL000218.1, GL000219.1, GL000220.1, GL000224.1, KI270438.1, KI270709.1, KI270719.1, KI270733.1, KI270738.1, KI270754.1, KI270757.1, KI270908.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, GL339449.2, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270728.1, KI270742.1, KI270745.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270871.1, KI270905.1, KZ208915.1, KZ559112.1, chrM, GL000218.1, GL000219.1, KI270719.1, KI270738.1, KI270754.1, KI270757.1, KI270908.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000225.1, KI270330.1, KI270442.1, KI270538.1, KI270580.1, KI270587.1, KI270724.1, KI270729.1, KI270736.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, GL339449.2, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270728.1, KI270742.1, KI270745.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270871.1, KI270905.1, KZ208915.1, KZ559112.1, chrM, GL000218.1, GL000219.1, KI270719.1, KI270738.1, KI270754.1, KI270757.1, KI270908.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000225.1, KI270330.1, KI270442.1, KI270538.1, KI270580.1, KI270587.1, KI270724.1, KI270729.1, KI270736.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270783.1, KI270784.1, KI270871.1, GL000218.1, GL000220.1, KI270709.1, KI270719.1, KI270733.1, KI270738.1, KI270757.1, KI270908.1, KN196487.1, GL000225.1, KI270330.1, KI270442.1, KI270538.1, KI270580.1, KI270587.1, KI270724.1
    ##   - in 'y': GL000253.2, GL000255.2, KI270712.1, KI270714.1, KI270721.1, KI270730.1, KI270761.1, KI270782.1, KI270832.1, KI270853.1, KI270857.1, KI270879.1, KI270880.1, KI270892.1, KI270897.1, KI270899.1, KN538364.1, KN538370.1, KN538372.1, KQ031384.1, KQ031389.1, KV575244.1, KZ208922.1, KZ559100.1, ML143345.1, ML143353.1, ML143355.1, ML143365.1, ML143375.1, ML143377.1, ML143379.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270783.1, KI270784.1, KI270871.1, GL000218.1, GL000220.1, KI270709.1, KI270719.1, KI270733.1, KI270738.1, KI270757.1, KI270908.1, KN196487.1, GL000225.1, KI270330.1, KI270442.1, KI270538.1, KI270580.1, KI270587.1, KI270724.1
    ##   - in 'y': GL000253.2, GL000255.2, KI270712.1, KI270714.1, KI270721.1, KI270730.1, KI270761.1, KI270782.1, KI270832.1, KI270853.1, KI270857.1, KI270879.1, KI270880.1, KI270892.1, KI270897.1, KI270899.1, KN538364.1, KN538370.1, KN538372.1, KQ031384.1, KQ031389.1, KV575244.1, KZ208922.1, KZ559100.1, ML143345.1, ML143353.1, ML143355.1, ML143365.1, ML143375.1, ML143377.1, ML143379.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KQ031384.1, ML143380.1
    ##   - in 'y': GL000220.1, GL000225.1, KI270320.1, KI270330.1, KI270538.1, KI270587.1, KI270709.1, KI270728.1, KI270744.1, chr10, chr14, chr15, chr16, chr18, chr22
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KQ031384.1, ML143380.1
    ##   - in 'y': GL000220.1, GL000225.1, KI270320.1, KI270330.1, KI270538.1, KI270587.1, KI270709.1, KI270728.1, KI270744.1, chr10, chr14, chr15, chr16, chr18, chr22
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000218.1, GL000221.1, GL000224.1, GL000225.1, GL000253.2, GL383563.3, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270711.1, KI270743.1, KI270745.1, KI270754.1, KI270762.1, KI270816.1, KI270830.1, KI270831.1, KI270853.1, KI270855.1, KI270856.1, KI270879.1, KI270905.1, KI270908.1, KV575244.1, KV766198.1, KV880768.1, KZ208915.1, ML143345.1, ML143355.1, ML143366.1
    ##   - in 'y': KI270337.1, KI270411.1, KI270466.1, KI270467.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000218.1, GL000221.1, GL000224.1, GL000225.1, GL000253.2, GL383563.3, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270711.1, KI270743.1, KI270745.1, KI270754.1, KI270762.1, KI270816.1, KI270830.1, KI270831.1, KI270853.1, KI270855.1, KI270856.1, KI270879.1, KI270905.1, KI270908.1, KV575244.1, KV766198.1, KV880768.1, KZ208915.1, ML143345.1, ML143355.1, ML143366.1
    ##   - in 'y': KI270337.1, KI270411.1, KI270466.1, KI270467.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL000256.2, KI270713.1, KI270714.1, KI270742.1, KI270905.1, KV766198.1, KZ208915.1, ML143372.1, ML143380.1, chr13, chr8
    ##   - in 'y': GL383563.3
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL000256.2, KI270713.1, KI270714.1, KI270742.1, KI270905.1, KV766198.1, KZ208915.1, ML143372.1, ML143380.1, chr13, chr8
    ##   - in 'y': GL383563.3
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL000256.2, KI270713.1, KI270714.1, KI270742.1, KI270905.1, KV766198.1, KZ208915.1, ML143372.1, ML143380.1, chr13, chr8
    ##   - in 'y': GL000216.2, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL000256.2, KI270713.1, KI270714.1, KI270742.1, KI270905.1, KV766198.1, KZ208915.1, ML143372.1, ML143380.1, chr13, chr8
    ##   - in 'y': GL000216.2, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270905.1, KZ208915.1, GL383563.3, GL000216.2
    ##   - in 'y': GL000214.1, GL000219.1, GL000253.2, GL949752.1, KI270442.1, KI270711.1, KI270718.1, KI270728.1, KI270733.1, KI270754.1, KI270761.1, KI270784.1, KI270816.1, KI270819.1, KI270850.1, KI270853.1, KI270857.1, KI270868.1, KI270869.1, KI270897.1, KI270899.1, KI270902.1, KI270908.1, KN538364.1, KN538370.1, KQ031389.1, KQ090026.1, KV575244.1, KV880764.1, KV880768.1, KZ559100.1, KZ559112.1, ML143345.1, ML143352.1, ML143353.1, ML143358.1, ML143365.1, ML143366.1, ML143371.1, ML143375.1, ML143377.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KI270905.1, KZ208915.1, GL383563.3, GL000216.2
    ##   - in 'y': GL000214.1, GL000219.1, GL000253.2, GL949752.1, KI270442.1, KI270711.1, KI270718.1, KI270728.1, KI270733.1, KI270754.1, KI270761.1, KI270784.1, KI270816.1, KI270819.1, KI270850.1, KI270853.1, KI270857.1, KI270868.1, KI270869.1, KI270897.1, KI270899.1, KI270902.1, KI270908.1, KN538364.1, KN538370.1, KQ031389.1, KQ090026.1, KV575244.1, KV880764.1, KV880768.1, KZ559100.1, KZ559112.1, ML143345.1, ML143352.1, ML143353.1, ML143358.1, ML143365.1, ML143366.1, ML143371.1, ML143375.1, ML143377.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000225.1, KI270466.1, KI270467.1, KI270538.1, KI270736.1, KI270765.1, KQ983257.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL000221.1, GL000251.2, GL000253.2, GL000256.2, GL339449.2, GL383522.1, GL383556.1, GL383580.2, GL949752.1, GL949753.2, JH159146.1, KI270707.1, KI270711.1, KI270712.1, KI270719.1, KI270723.1, KI270726.1, KI270731.1, KI270733.1, KI270734.1, KI270743.1, KI270745.1, KI270751.1, KI270762.1, KI270770.1, KI270780.1, KI270782.1, KI270783.1, KI270784.1, KI270787.1, KI270792.1, KI270803.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270832.1, KI270836.1, KI270842.1, KI270850.1, KI270857.1, KI270860.1, KI270861.1, KI270876.1, KI270880.1, KI270902.1, KI270903.1, KN196472.1, KN196484.1, KN538361.1, KN538372.1, KQ458383.1, KQ458385.1, KV766193.1, KV766196.1, KV880763.1, KV880764.1, KZ208912.1, KZ208917.1, KZ208921.1, KZ559103.1, KZ559105.1, KZ559112.1, ML143355.1, ML143367.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000225.1, KI270466.1, KI270467.1, KI270538.1, KI270736.1, KI270765.1, KQ983257.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL000221.1, GL000251.2, GL000253.2, GL000256.2, GL339449.2, GL383522.1, GL383556.1, GL383580.2, GL949752.1, GL949753.2, JH159146.1, KI270707.1, KI270711.1, KI270712.1, KI270719.1, KI270723.1, KI270726.1, KI270731.1, KI270733.1, KI270734.1, KI270743.1, KI270745.1, KI270751.1, KI270762.1, KI270770.1, KI270780.1, KI270782.1, KI270783.1, KI270784.1, KI270787.1, KI270792.1, KI270803.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270832.1, KI270836.1, KI270842.1, KI270850.1, KI270857.1, KI270860.1, KI270861.1, KI270876.1, KI270880.1, KI270902.1, KI270903.1, KN196472.1, KN196484.1, KN538361.1, KN538372.1, KQ458383.1, KQ458385.1, KV766193.1, KV766196.1, KV880763.1, KV880764.1, KZ208912.1, KZ208917.1, KZ208921.1, KZ559103.1, KZ559105.1, KZ559112.1, ML143355.1, ML143367.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270330.1, KI270712.1, KI270735.1, KI270751.1, KI270782.1, KQ458384.1, KV880764.1, KV880768.1, KZ208913.1, ML143344.1, ML143345.1, ML143352.1
    ##   - in 'y': GL000218.1, GL000224.1, GL000254.2, GL383527.1, KI270709.1, KI270713.1, KI270714.1, KI270745.1, KI270805.1, KI270831.1, KI270857.1, KI270904.1, KN196487.1, KN538361.1, KN538372.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270330.1, KI270712.1, KI270735.1, KI270751.1, KI270782.1, KQ458384.1, KV880764.1, KV880768.1, KZ208913.1, ML143344.1, ML143345.1, ML143352.1
    ##   - in 'y': GL000218.1, GL000224.1, GL000254.2, GL383527.1, KI270709.1, KI270713.1, KI270714.1, KI270745.1, KI270805.1, KI270831.1, KI270857.1, KI270904.1, KN196487.1, KN538361.1, KN538372.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270706.1, KI270742.1, KI270765.1, KZ208915.1, ML143345.1, ML143377.1
    ##   - in 'y': GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270706.1, KI270742.1, KI270765.1, KZ208915.1, ML143345.1, ML143377.1
    ##   - in 'y': GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL000256.2, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270765.1, KI270850.1, KI270879.1, KI270908.1, KV575244.1, KV766198.1, KZ208915.1, ML143345.1, ML143377.1, GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   - in 'y': GL000009.2, KI270442.1, KI270787.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL000256.2, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270765.1, KI270850.1, KI270879.1, KI270908.1, KV575244.1, KV766198.1, KZ208915.1, ML143345.1, ML143377.1, GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   - in 'y': GL000009.2, KI270442.1, KI270787.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL000256.2, KI270711.1, KI270714.1, KI270765.1, KI270850.1, KI270879.1, KI270908.1, KV575244.1, KV766198.1, ML143345.1, ML143377.1, GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM, GL000009.2, KI270787.1
    ##   - in 'y': GL000219.1, GL877875.1, KI270878.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL000256.2, KI270711.1, KI270714.1, KI270765.1, KI270850.1, KI270879.1, KI270908.1, KV575244.1, KV766198.1, ML143345.1, ML143377.1, GL000251.2, KI270905.1, KQ458383.1, KV880768.1, ML143375.1, chrM, GL000009.2, KI270787.1
    ##   - in 'y': GL000219.1, GL877875.1, KI270878.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000214.1, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL339449.2, GL383520.2, GL383522.1, GL877875.1, GL949746.1, JH159146.1, KI270707.1, KI270721.1, KI270728.1, KI270734.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270792.1, KI270803.1, KI270819.1, KI270821.1, KI270831.1, KI270848.1, KI270850.1, KI270851.1, KI270856.1, KI270857.1, KI270861.1, KI270868.1, KI270878.1, KI270880.1, KI270904.1, KI270905.1, KN538361.1, KN538364.1, KQ031389.1, KQ458385.1, KV880764.1, KV880768.1, KZ208907.1, KZ208912.1, KZ208913.1, KZ208921.1, KZ559105.1, ML143345.1, ML143353.1, ML143355.1, ML143365.1, ML143366.1, ML143372.1, ML143375.1, ML143379.1, ML143380.1
    ##   - in 'y': GL949752.1, KI270767.1, KI270782.1, KI270893.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000214.1, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL339449.2, GL383520.2, GL383522.1, GL877875.1, GL949746.1, JH159146.1, KI270707.1, KI270721.1, KI270728.1, KI270734.1, KI270765.1, KI270783.1, KI270784.1, KI270787.1, KI270792.1, KI270803.1, KI270819.1, KI270821.1, KI270831.1, KI270848.1, KI270850.1, KI270851.1, KI270856.1, KI270857.1, KI270861.1, KI270868.1, KI270878.1, KI270880.1, KI270904.1, KI270905.1, KN538361.1, KN538364.1, KQ031389.1, KQ458385.1, KV880764.1, KV880768.1, KZ208907.1, KZ208912.1, KZ208913.1, KZ208921.1, KZ559105.1, ML143345.1, ML143353.1, ML143355.1, ML143365.1, ML143366.1, ML143372.1, ML143375.1, ML143379.1, ML143380.1
    ##   - in 'y': GL949752.1, KI270767.1, KI270782.1, KI270893.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000254.2, JH159146.1, KI270792.1, KI270821.1, KI270904.1, KZ208913.1, KI270767.1, KI270893.1
    ##   - in 'y': GL383567.1, GL383581.2, KI270717.1, KI270726.1, KI270754.1, KI270876.1, KI270897.1, KI270899.1, KI270934.1, KI270936.1, KI270937.1, KN538360.1, KN538370.1, KQ090016.1, KZ208922.1, ML143350.1, ML143352.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000254.2, JH159146.1, KI270792.1, KI270821.1, KI270904.1, KZ208913.1, KI270767.1, KI270893.1
    ##   - in 'y': GL383567.1, GL383581.2, KI270717.1, KI270726.1, KI270754.1, KI270876.1, KI270897.1, KI270899.1, KI270934.1, KI270936.1, KI270937.1, KN538360.1, KN538370.1, KQ090016.1, KZ208922.1, ML143350.1, ML143352.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000214.1, GL000219.1, GL877875.1, JH159146.1, KI270803.1, KI270819.1, KI270821.1, KI270878.1, KI270880.1, KZ208907.1, KZ208913.1, KZ208921.1, ML143365.1, ML143366.1, GL383567.1, GL383581.2, KI270717.1, KI270726.1, KI270876.1, KI270899.1, KI270934.1, KI270936.1, KI270937.1, KN538360.1, KZ208922.1, ML143350.1, ML143352.1, ML143358.1
    ##   - in 'y': GL383578.2, KI270731.1, KI270760.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000214.1, GL000219.1, GL877875.1, JH159146.1, KI270803.1, KI270819.1, KI270821.1, KI270878.1, KI270880.1, KZ208907.1, KZ208913.1, KZ208921.1, ML143365.1, ML143366.1, GL383567.1, GL383581.2, KI270717.1, KI270726.1, KI270876.1, KI270899.1, KI270934.1, KI270936.1, KI270937.1, KN538360.1, KZ208922.1, ML143350.1, ML143352.1, ML143358.1
    ##   - in 'y': GL383578.2, KI270731.1, KI270760.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, KI270538.1
    ##   - in 'y': GL000008.2, GL000219.1, GL000221.1, GL000253.2, GL000254.2, GL000256.2, GL949752.1, KI270337.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270719.1, KI270741.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270770.1, KI270784.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270836.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270900.1, KI270905.1, KI270908.1, KI270936.1, KI270937.1, KQ090026.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559109.1, ML143366.1, ML143367.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, KI270538.1
    ##   - in 'y': GL000008.2, GL000219.1, GL000221.1, GL000253.2, GL000254.2, GL000256.2, GL949752.1, KI270337.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270719.1, KI270741.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270770.1, KI270784.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270836.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270900.1, KI270905.1, KI270908.1, KI270936.1, KI270937.1, KQ090026.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559109.1, ML143366.1, ML143367.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, KI270538.1, GL000008.2, GL000221.1, GL000254.2, KI270337.1, KI270466.1, KI270467.1, KI270711.1, KI270827.1, KI270830.1, KI270856.1, KI270880.1, KI270937.1, KV880768.1, KZ208912.1, KZ559109.1, chrM
    ##   - in 'y': GL000214.1, KI270765.1, KI270897.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, KI270538.1, GL000008.2, GL000221.1, GL000254.2, KI270337.1, KI270466.1, KI270467.1, KI270711.1, KI270827.1, KI270830.1, KI270856.1, KI270880.1, KI270937.1, KV880768.1, KZ208912.1, KZ559109.1, chrM
    ##   - in 'y': GL000214.1, KI270765.1, KI270897.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000221.1, GL000253.2, GL000254.2, GL000256.2, GL949752.1, KI270337.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270719.1, KI270741.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270770.1, KI270784.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270836.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270900.1, KI270905.1, KI270908.1, KI270936.1, KI270937.1, KQ090026.1, KV880768.1, KZ208912.1, KZ559109.1, ML143366.1, ML143367.1, ML143372.1, chrM, GL000214.1, KI270765.1, KI270897.1
    ##   - in 'y': KI270438.1, KI270709.1, KQ031389.1, KZ208907.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000221.1, GL000253.2, GL000254.2, GL000256.2, GL949752.1, KI270337.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270719.1, KI270741.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270770.1, KI270784.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270836.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270900.1, KI270905.1, KI270908.1, KI270936.1, KI270937.1, KQ090026.1, KV880768.1, KZ208912.1, KZ559109.1, ML143366.1, ML143367.1, ML143372.1, chrM, GL000214.1, KI270765.1, KI270897.1
    ##   - in 'y': KI270438.1, KI270709.1, KQ031389.1, KZ208907.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, GL949752.1, KI270765.1, KI270857.1, KI270905.1, KN196479.1, KV880768.1, ML143366.1
    ##   - in 'y': GL000219.1, GL339449.2, KI270466.1, KI270467.1, KI270706.1, KI270714.1, KI270744.1, KI270822.1, KI270832.1, KI270849.1, KI270880.1, KI270908.1, KQ458383.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, GL949752.1, KI270765.1, KI270857.1, KI270905.1, KN196479.1, KV880768.1, ML143366.1
    ##   - in 'y': GL000219.1, GL339449.2, KI270466.1, KI270467.1, KI270706.1, KI270714.1, KI270744.1, KI270822.1, KI270832.1, KI270849.1, KI270880.1, KI270908.1, KQ458383.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000252.2, GL000254.2, GL339449.2, GL383574.1, GL383578.2, GL383579.2, GL877875.1, GL949746.1, KI270442.1, KI270707.1, KI270721.1, KI270728.1, KI270731.1, KI270745.1, KI270765.1, KI270783.1, KI270787.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270857.1, KI270862.1, KI270904.1, KI270908.1, KN538364.1, KQ031389.1, KQ458383.1, KQ458384.1, KV880764.1, KZ208912.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143345.1, ML143352.1
    ##   - in 'y': chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000252.2, GL000254.2, GL339449.2, GL383574.1, GL383578.2, GL383579.2, GL877875.1, GL949746.1, KI270442.1, KI270707.1, KI270721.1, KI270728.1, KI270731.1, KI270745.1, KI270765.1, KI270783.1, KI270787.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270857.1, KI270862.1, KI270904.1, KI270908.1, KN538364.1, KQ031389.1, KQ458383.1, KQ458384.1, KV880764.1, KZ208912.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143345.1, ML143352.1
    ##   - in 'y': chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, GL383563.3, GL383580.2, GL383581.2, KI270438.1, KI270723.1, KI270804.1, KI270936.1, KN538372.1, KZ559105.1, KZ559109.1, ML143366.1
    ##   - in 'y': GL000251.2, KI270731.1, KI270734.1, KI270737.1, KI270749.1, KI270831.1, KI270837.1, KI270848.1, KN538361.1, KQ458383.1, KV880764.1, ML143344.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, GL383563.3, GL383580.2, GL383581.2, KI270438.1, KI270723.1, KI270804.1, KI270936.1, KN538372.1, KZ559105.1, KZ559109.1, ML143366.1
    ##   - in 'y': GL000251.2, KI270731.1, KI270734.1, KI270737.1, KI270749.1, KI270831.1, KI270837.1, KI270848.1, KN538361.1, KQ458383.1, KV880764.1, ML143344.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000253.2, GL339449.2, GL383556.1, KI270442.1, KI270720.1, KI270737.1, KI270750.1, KI270802.1, KI270821.1, KI270822.1, KI270853.1, KI270896.1, KV766198.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000251.2, KI270746.1, KN196487.1, KN538372.1, KQ983257.1, KV880768.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000253.2, GL339449.2, GL383556.1, KI270442.1, KI270720.1, KI270737.1, KI270750.1, KI270802.1, KI270821.1, KI270822.1, KI270853.1, KI270896.1, KV766198.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000251.2, KI270746.1, KN196487.1, KN538372.1, KQ983257.1, KV880768.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270311.1, KI270315.1, KI270435.1, KI270508.1, KI270591.1, KI270708.1, KI270720.1, KI270813.1, KI270831.1, KI270879.1, KN196472.1, KN538372.1, KQ031384.1, KQ458383.1, KV766192.1, KZ208919.1
    ##   - in 'y': GL000253.2, GL383526.1, GL383542.1, GL383545.1, KI270519.1, KI270731.1, KI270734.1, KI270821.1, KI270830.1, KI270899.1, KI270904.1, KI270934.1, KN196480.1, KQ458385.1, KZ208912.1, KZ208921.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270311.1, KI270315.1, KI270435.1, KI270508.1, KI270591.1, KI270708.1, KI270720.1, KI270813.1, KI270831.1, KI270879.1, KN196472.1, KN538372.1, KQ031384.1, KQ458383.1, KV766192.1, KZ208919.1
    ##   - in 'y': GL000253.2, GL383526.1, GL383542.1, GL383545.1, KI270519.1, KI270731.1, KI270734.1, KI270821.1, KI270830.1, KI270899.1, KI270904.1, KI270934.1, KN196480.1, KQ458385.1, KZ208912.1, KZ208921.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270589.1, KI270747.1, KI270761.1, KI270810.1, KI270879.1, KI270892.1
    ##   - in 'y': GL000252.2, GL339449.2, GL383527.1, GL383542.1, JH159136.1, KI270442.1, KI270706.1, KI270707.1, KI270714.1, KI270722.1, KI270746.1, KI270783.1, KI270802.1, KI270804.1, KI270809.1, KI270830.1, KI270836.1, KI270851.1, KN196484.1, KQ031389.1, KV766192.1, KV880763.1, ML143378.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270589.1, KI270747.1, KI270761.1, KI270810.1, KI270879.1, KI270892.1
    ##   - in 'y': GL000252.2, GL339449.2, GL383527.1, GL383542.1, JH159136.1, KI270442.1, KI270706.1, KI270707.1, KI270714.1, KI270722.1, KI270746.1, KI270783.1, KI270802.1, KI270804.1, KI270809.1, KI270830.1, KI270836.1, KI270851.1, KN196484.1, KQ031389.1, KV766192.1, KV880763.1, ML143378.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000255.2, KI270538.1, KI270745.1, KI270754.1, KI270787.1, KI270803.1, ML143367.1
    ##   - in 'y': GL383567.1, GL949752.1, KI270411.1, KI270438.1, KI270728.1, KI270819.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000255.2, KI270538.1, KI270745.1, KI270754.1, KI270787.1, KI270803.1, ML143367.1
    ##   - in 'y': GL383567.1, GL949752.1, KI270411.1, KI270438.1, KI270728.1, KI270819.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, KI270787.1
    ##   - in 'y': GL000219.1, GL877875.1, KI270706.1, KI270713.1, KI270878.1, KZ208915.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, KI270787.1
    ##   - in 'y': GL000219.1, GL877875.1, KI270706.1, KI270713.1, KI270878.1, KZ208915.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270706.1, KI270754.1, KI270849.1
    ##   - in 'y': GL000194.1, GL000221.1, GL383563.3, GL383567.1, KI270714.1, KI270724.1, KI270728.1, KI270744.1, KI270830.1, KI270897.1, KI270924.1, KN538372.1, KQ983257.1, KV880768.1, ML143353.1, ML143367.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270706.1, KI270754.1, KI270849.1
    ##   - in 'y': GL000194.1, GL000221.1, GL383563.3, GL383567.1, KI270714.1, KI270724.1, KI270728.1, KI270744.1, KI270830.1, KI270897.1, KI270924.1, KN538372.1, KQ983257.1, KV880768.1, ML143353.1, ML143367.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1
    ##   - in 'y': KI270337.1, KI270713.1, chr10, chr18, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1
    ##   - in 'y': KI270337.1, KI270713.1, chr10, chr18, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000225.1, GL000251.2, GL000254.2, GL339449.2, KI270538.1, KI270734.1, KI270821.1, KI270908.1
    ##   - in 'y': GL000214.1, KI270779.1, KI270782.1, KI270787.1, KI270830.1, KI270832.1, KI270848.1, KI270850.1, KI270853.1, KI270924.1, KV575244.1, KV880763.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000225.1, GL000251.2, GL000254.2, GL339449.2, KI270538.1, KI270734.1, KI270821.1, KI270908.1
    ##   - in 'y': GL000214.1, KI270779.1, KI270782.1, KI270787.1, KI270830.1, KI270832.1, KI270848.1, KI270850.1, KI270853.1, KI270924.1, KV575244.1, KV880763.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270544.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000224.1, GL383522.1, KI270304.1, KI270732.1, KI270744.1, KI270754.1, KI270878.1, KQ031384.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270544.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000224.1, GL383522.1, KI270304.1, KI270732.1, KI270744.1, KI270754.1, KI270878.1, KQ031384.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270782.1, KI270819.1, KQ090026.1
    ##   - in 'y': GL000255.2, KI270728.1, KI270731.1, KI270783.1, KI270816.1, KV575244.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270782.1, KI270819.1, KQ090026.1
    ##   - in 'y': GL000255.2, KI270728.1, KI270731.1, KI270783.1, KI270816.1, KV575244.1, ML143367.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, KI270713.1, KI270737.1, KI270751.1, KZ208915.1, ML143377.1, ML143380.1
    ##   - in 'y': GL000225.1, GL339449.2, KI270764.1, KI270908.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, KI270713.1, KI270737.1, KI270751.1, KZ208915.1, ML143377.1, ML143380.1
    ##   - in 'y': GL000225.1, GL339449.2, KI270764.1, KI270908.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000208.1, GL000221.1, GL339449.2, KI270320.1, KI270411.1, KI270538.1, KI270544.1, KI270723.1, KI270849.1, KI270892.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383574.1, GL383578.2, KI270465.1, KI270730.1, KI270731.1, KI270899.1, KN196484.1, KN538360.1, KQ031384.1, KQ983257.1, KV575244.1, KV880764.1, ML143360.1, ML143366.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000208.1, GL000221.1, GL339449.2, KI270320.1, KI270411.1, KI270538.1, KI270544.1, KI270723.1, KI270849.1, KI270892.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383574.1, GL383578.2, KI270465.1, KI270730.1, KI270731.1, KI270899.1, KN196484.1, KN538360.1, KQ031384.1, KQ983257.1, KV575244.1, KV880764.1, ML143360.1, ML143366.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, KN196487.1
    ##   - in 'y': GL000253.2, KI270310.1, KI270438.1, KI270519.1, KI270544.1, KI270706.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, KN196487.1
    ##   - in 'y': GL000253.2, KI270310.1, KI270438.1, KI270519.1, KI270544.1, KI270706.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000224.1, KI270304.1, KI270311.1, KI270320.1, KI270737.1
    ##   - in 'y': GL000205.2, GL000251.2, GL383555.2, GL383567.1, GL877875.1, KI270711.1, KI270713.1, KI270728.1, KI270783.1, KI270787.1, KI270803.1, KI270821.1, KI270878.1, KI270908.1, KQ031389.1, KQ090016.1, KV766198.1, KV880768.1, KZ559112.1, ML143367.1, ML143371.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000224.1, KI270304.1, KI270311.1, KI270320.1, KI270737.1
    ##   - in 'y': GL000205.2, GL000251.2, GL383555.2, GL383567.1, GL877875.1, KI270711.1, KI270713.1, KI270728.1, KI270783.1, KI270787.1, KI270803.1, KI270821.1, KI270878.1, KI270908.1, KQ031389.1, KQ090016.1, KV766198.1, KV880768.1, KZ559112.1, ML143367.1, ML143371.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000218.1, GL000255.2, KI270733.1, KI270754.1, KI270853.1, KI270892.1, ML143372.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000256.2, JH159146.1, KI270712.1, KI270772.1, KQ983257.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000218.1, GL000255.2, KI270733.1, KI270754.1, KI270853.1, KI270892.1, ML143372.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000256.2, JH159146.1, KI270712.1, KI270772.1, KQ983257.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000255.2, GL339449.2, KI270707.1, KI270712.1, KI270714.1, KI270782.1, KI270830.1, KI270836.1, KI270856.1, KI270895.1, KI270908.1, KN196484.1, KN196487.1, KQ031384.1, KQ983257.1, KV880768.1, KZ208915.1, ML143372.1, ML143377.1
    ##   - in 'y': KI270442.1, KI270723.1, KI270804.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000255.2, GL339449.2, KI270707.1, KI270712.1, KI270714.1, KI270782.1, KI270830.1, KI270836.1, KI270856.1, KI270895.1, KI270908.1, KN196484.1, KN196487.1, KQ031384.1, KQ983257.1, KV880768.1, KZ208915.1, ML143372.1, ML143377.1
    ##   - in 'y': KI270442.1, KI270723.1, KI270804.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000216.2, KI270707.1, KI270712.1, KI270728.1, KI270731.1, KI270836.1, KI270895.1, KI270908.1, KN196484.1, KN196487.1, KQ031384.1, KQ983257.1, KI270723.1, KI270804.1
    ##   - in 'y': GL000256.2, KI270310.1, KI270706.1, KI270736.1, KI270743.1, KI270761.1, KI270770.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KI270905.1, KQ090026.1, KV575244.1, KV766198.1, ML143345.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000216.2, KI270707.1, KI270712.1, KI270728.1, KI270731.1, KI270836.1, KI270895.1, KI270908.1, KN196484.1, KN196487.1, KQ031384.1, KQ983257.1, KI270723.1, KI270804.1
    ##   - in 'y': GL000256.2, KI270310.1, KI270706.1, KI270736.1, KI270743.1, KI270761.1, KI270770.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KI270905.1, KQ090026.1, KV575244.1, KV766198.1, ML143345.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000216.2, GL000225.1, GL339449.2, KI270707.1, KI270712.1, KI270728.1, KI270782.1, KI270836.1, KI270879.1, KI270895.1, KI270908.1, KN196487.1, KQ031384.1, KQ983257.1, KZ208915.1, ML143377.1, KI270723.1, KI270804.1, KI270706.1, KI270736.1, KI270761.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KQ090026.1, KV575244.1
    ##   - in 'y': ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000216.2, GL000225.1, GL339449.2, KI270707.1, KI270712.1, KI270728.1, KI270782.1, KI270836.1, KI270879.1, KI270895.1, KI270908.1, KN196487.1, KQ031384.1, KQ983257.1, KZ208915.1, ML143377.1, KI270723.1, KI270804.1, KI270706.1, KI270736.1, KI270761.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KQ090026.1, KV575244.1
    ##   - in 'y': ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000225.1, GL339449.2, KI270706.1, KI270736.1, KI270761.1, KI270782.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KI270879.1, KQ090026.1, KV575244.1, KZ208915.1, ML143377.1
    ##   - in 'y': KI270731.1, KN196484.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000225.1, GL339449.2, KI270706.1, KI270736.1, KI270761.1, KI270782.1, KI270783.1, KI270819.1, KI270831.1, KI270868.1, KI270879.1, KQ090026.1, KV575244.1, KZ208915.1, ML143377.1
    ##   - in 'y': KI270731.1, KN196484.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, KI270711.1, KI270734.1, KI270750.1, KI270754.1, KI270765.1, KI270787.1, KI270818.1, KI270856.1, KI270857.1, KQ090026.1
    ##   - in 'y': GL000218.1, GL000220.1, GL000224.1, GL339449.2, KI270310.1, KI270708.1, KI270712.1, KI270714.1, KI270722.1, KI270733.1, KI270743.1, KI270746.1, KI270763.1, KI270772.1, KI270832.1, KI270878.1, KN538361.1, KN538364.1, KQ031389.1, KQ090027.1, KQ458383.1, KV575244.1, KZ208922.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, KI270711.1, KI270734.1, KI270750.1, KI270754.1, KI270765.1, KI270787.1, KI270818.1, KI270856.1, KI270857.1, KQ090026.1
    ##   - in 'y': GL000218.1, GL000220.1, GL000224.1, GL339449.2, KI270310.1, KI270708.1, KI270712.1, KI270714.1, KI270722.1, KI270733.1, KI270743.1, KI270746.1, KI270763.1, KI270772.1, KI270832.1, KI270878.1, KN538361.1, KN538364.1, KQ031389.1, KQ090027.1, KQ458383.1, KV575244.1, KZ208922.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383542.1, KI270712.1, KI270784.1, KQ090026.1, KV766198.1
    ##   - in 'y': GL000218.1, GL000224.1, GL383567.1, GL877875.1, KI270733.1, KI270743.1, KI270782.1, KI270819.1, KQ031389.1, KV880765.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383542.1, KI270712.1, KI270784.1, KQ090026.1, KV766198.1
    ##   - in 'y': GL000218.1, GL000224.1, GL383567.1, GL877875.1, KI270733.1, KI270743.1, KI270782.1, KI270819.1, KQ031389.1, KV880765.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, KI270787.1, KN538364.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270731.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, KI270787.1, KN538364.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270731.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000205.2, GL000218.1, GL000221.1, GL000224.1, GL000254.2, GL339449.2, GL383527.1, GL383581.2, JH159146.1, KI270707.1, KI270720.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270734.1, KI270743.1, KI270750.1, KI270751.1, KI270765.1, KI270772.1, KI270783.1, KI270787.1, KI270803.1, KI270810.1, KI270811.1, KI270821.1, KI270822.1, KI270856.1, KI270857.1, KI270869.1, KI270878.1, KI270896.1, KI270905.1, KN538364.1, KQ090016.1, KQ090026.1, KQ983255.1, KV575244.1, KZ208913.1, ML143372.1, ML143377.1, GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270731.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1
    ##   - in 'y': KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270733.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000205.2, GL000218.1, GL000221.1, GL000224.1, GL000254.2, GL339449.2, GL383527.1, GL383581.2, JH159146.1, KI270707.1, KI270720.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270734.1, KI270743.1, KI270750.1, KI270751.1, KI270765.1, KI270772.1, KI270783.1, KI270787.1, KI270803.1, KI270810.1, KI270811.1, KI270821.1, KI270822.1, KI270856.1, KI270857.1, KI270869.1, KI270878.1, KI270896.1, KI270905.1, KN538364.1, KQ090016.1, KQ090026.1, KQ983255.1, KV575244.1, KZ208913.1, ML143372.1, ML143377.1, GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270731.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1
    ##   - in 'y': KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270733.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000194.1, GL000205.2, GL000218.1, GL000219.1, GL000221.1, GL000254.2, GL339449.2, GL383527.1, GL383581.2, JH159146.1, KI270706.1, KI270707.1, KI270713.1, KI270720.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270734.1, KI270743.1, KI270744.1, KI270750.1, KI270751.1, KI270765.1, KI270772.1, KI270783.1, KI270787.1, KI270803.1, KI270810.1, KI270811.1, KI270821.1, KI270822.1, KI270849.1, KI270850.1, KI270856.1, KI270857.1, KI270869.1, KI270878.1, KI270879.1, KI270896.1, KI270905.1, KN538364.1, KQ090016.1, KQ090026.1, KQ983255.1, KV575244.1, KV880768.1, KZ208913.1, ML143355.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1, chrM, GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1, KI270337.1, KI270466.1, KI270467.1, KZ559112.1
    ##   - in 'y': KI270310.1, KI270538.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000194.1, GL000205.2, GL000218.1, GL000219.1, GL000221.1, GL000254.2, GL339449.2, GL383527.1, GL383581.2, JH159146.1, KI270706.1, KI270707.1, KI270713.1, KI270720.1, KI270723.1, KI270724.1, KI270725.1, KI270728.1, KI270734.1, KI270743.1, KI270744.1, KI270750.1, KI270751.1, KI270765.1, KI270772.1, KI270783.1, KI270787.1, KI270803.1, KI270810.1, KI270811.1, KI270821.1, KI270822.1, KI270849.1, KI270850.1, KI270856.1, KI270857.1, KI270869.1, KI270878.1, KI270879.1, KI270896.1, KI270905.1, KN538364.1, KQ090016.1, KQ090026.1, KQ983255.1, KV575244.1, KV880768.1, KZ208913.1, ML143355.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1, chrM, GL000214.1, GL000225.1, GL000257.2, KI270708.1, KI270714.1, KI270729.1, KI270735.1, KI270745.1, KI270748.1, KI270760.1, KI270808.1, KI270827.1, KI270830.1, KI270926.1, KV766198.1, KZ559103.1, ML143364.1, KI270337.1, KI270466.1, KI270467.1, KZ559112.1
    ##   - in 'y': KI270310.1, KI270538.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000219.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270744.1, KI270849.1, KI270850.1, KI270879.1, KV880768.1, KZ559112.1, ML143355.1, ML143366.1, ML143380.1, chrM
    ##   - in 'y': GL000224.1, KI270310.1, KI270538.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000219.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270744.1, KI270849.1, KI270850.1, KI270879.1, KV880768.1, KZ559112.1, ML143355.1, ML143366.1, ML143380.1, chrM
    ##   - in 'y': GL000224.1, KI270310.1, KI270538.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000251.2, KI270467.1, KI270719.1, KI270721.1, KI270878.1, ML143369.1
    ##   - in 'y': GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270830.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000251.2, KI270467.1, KI270719.1, KI270721.1, KI270878.1, ML143369.1
    ##   - in 'y': GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270830.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, KI270711.1, KI270731.1, KI270810.1, KI270853.1, KI270876.1
    ##   - in 'y': GL000221.1, GL000224.1, KI270743.1, KI270767.1, KI270782.1, KI270787.1, KI270816.1, KI270830.1, KI270831.1, KI270880.1, KI270892.1, KV880764.1, KZ559105.1, ML143345.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, KI270711.1, KI270731.1, KI270810.1, KI270853.1, KI270876.1
    ##   - in 'y': GL000221.1, GL000224.1, KI270743.1, KI270767.1, KI270782.1, KI270787.1, KI270816.1, KI270830.1, KI270831.1, KI270880.1, KI270892.1, KV880764.1, KZ559105.1, ML143345.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, KI270733.1, KI270821.1, KQ090026.1
    ##   - in 'y': GL000253.2, GL339449.2, GL383522.1, KI270734.1, KI270827.1, KI270849.1, KI270873.1, KI270913.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, KI270733.1, KI270821.1, KQ090026.1
    ##   - in 'y': GL000253.2, GL339449.2, GL383522.1, KI270734.1, KI270827.1, KI270849.1, KI270873.1, KI270913.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270336.1, KI270908.1
    ##   - in 'y': GL000009.2, GL000224.1, GL339449.2, GL383556.1, GL383581.2, KI270515.1, KI270591.1, KI270706.1, KI270712.1, KI270720.1, KI270732.1, KI270743.1, KI270745.1, KI270746.1, KI270784.1, KI270821.1, KI270832.1, KI270849.1, KI270857.1, KI270904.1, KN538372.1, KN538373.1, KQ458383.1, KZ208907.1, KZ208912.1, KZ208913.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270336.1, KI270908.1
    ##   - in 'y': GL000009.2, GL000224.1, GL339449.2, GL383556.1, GL383581.2, KI270515.1, KI270591.1, KI270706.1, KI270712.1, KI270720.1, KI270732.1, KI270743.1, KI270745.1, KI270746.1, KI270784.1, KI270821.1, KI270832.1, KI270849.1, KI270857.1, KI270904.1, KN538372.1, KN538373.1, KQ458383.1, KZ208907.1, KZ208912.1, KZ208913.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270745.1, ML143371.1, ML143372.1
    ##   - in 'y': GL000218.1, GL000220.1, KI270709.1, KI270720.1, KQ458383.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270745.1, ML143371.1, ML143372.1
    ##   - in 'y': GL000218.1, GL000220.1, KI270709.1, KI270720.1, KQ458383.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000205.2, GL877875.1, KI270438.1
    ##   - in 'y': GL000008.2, GL000250.2, KI270713.1, KI270745.1, KI270783.1, KI270905.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000205.2, GL877875.1, KI270438.1
    ##   - in 'y': GL000008.2, GL000250.2, KI270713.1, KI270745.1, KI270783.1, KI270905.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1
    ##   - in 'y': GL000218.1, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1
    ##   - in 'y': GL000218.1, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270857.1, KQ458383.1, ML143372.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000255.2, GL877875.1, KI270538.1, KI270754.1, KI270805.1, KI270830.1, KI270853.1, KQ090016.1, KZ208915.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270857.1, KQ458383.1, ML143372.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000255.2, GL877875.1, KI270538.1, KI270754.1, KI270805.1, KI270830.1, KI270853.1, KQ090016.1, KZ208915.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383567.1, KI270706.1, KI270709.1, KI270724.1, KI270729.1, KI270754.1, KI270905.1, KQ031384.1, KV766198.1, ML143367.1
    ##   - in 'y': GL000224.1, KI270712.1, KI270733.1, KI270764.1, KI270782.1, KI270819.1, KI270830.1, KI270850.1, KN538364.1, ML143355.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383567.1, KI270706.1, KI270709.1, KI270724.1, KI270729.1, KI270754.1, KI270905.1, KQ031384.1, KV766198.1, ML143367.1
    ##   - in 'y': GL000224.1, KI270712.1, KI270733.1, KI270764.1, KI270782.1, KI270819.1, KI270830.1, KI270850.1, KN538364.1, ML143355.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000220.1, KI270337.1, KI270729.1, KI270733.1, KI270803.1
    ##   - in 'y': GL000219.1, KI270713.1, KI270844.1, KI270908.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000220.1, KI270337.1, KI270729.1, KI270733.1, KI270803.1
    ##   - in 'y': GL000219.1, KI270713.1, KI270844.1, KI270908.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270442.1, KI270744.1, KI270782.1, KN196479.1
    ##   - in 'y': KI270764.1, chr14, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270442.1, KI270744.1, KI270782.1, KN196479.1
    ##   - in 'y': KI270764.1, chr14, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270803.1, KI270819.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270803.1, KI270819.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270322.1, KI270509.1, KI270517.1, KI270724.1, KI270735.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000208.1, GL000221.1, GL000256.2, GL383550.2, KI270304.1, KI270320.1, KI270515.1, KI270583.1, KI270584.1, KI270587.1, KI270593.1, KI270707.1, KI270729.1, KI270743.1, KI270772.1, KI270822.1, KI270908.1, KQ759759.1, KZ559103.1, KZ559110.1, ML143356.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270322.1, KI270509.1, KI270517.1, KI270724.1, KI270735.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000208.1, GL000221.1, GL000256.2, GL383550.2, KI270304.1, KI270320.1, KI270515.1, KI270583.1, KI270584.1, KI270587.1, KI270593.1, KI270707.1, KI270729.1, KI270743.1, KI270772.1, KI270822.1, KI270908.1, KQ759759.1, KZ559103.1, KZ559110.1, ML143356.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383545.1, GL383556.1, KI270333.1, KI270336.1, KI270544.1, KI270711.1, KI270714.1, KI270719.1, KI270731.1, KI270806.1, KI270818.1, KI270821.1, KI270925.1, KZ559105.1, ML143371.1
    ##   - in 'y': GL000216.2, GL383531.1, GL383578.2, KI270330.1, KI270411.1, KI270722.1, KI270741.1, KI270857.1, KI270871.1, KN196487.1, KQ090026.1, KZ559112.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383545.1, GL383556.1, KI270333.1, KI270336.1, KI270544.1, KI270711.1, KI270714.1, KI270719.1, KI270731.1, KI270806.1, KI270818.1, KI270821.1, KI270925.1, KZ559105.1, ML143371.1
    ##   - in 'y': GL000216.2, GL383531.1, GL383578.2, KI270330.1, KI270411.1, KI270722.1, KI270741.1, KI270857.1, KI270871.1, KN196487.1, KQ090026.1, KZ559112.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000214.1, GL000220.1, GL000255.2, KI270713.1, KI270743.1, KI270745.1, KV575244.1, ML143372.1, ML143378.1
    ##   - in 'y': GL000225.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000214.1, GL000220.1, GL000255.2, KI270713.1, KI270743.1, KI270745.1, KV575244.1, ML143372.1, ML143378.1
    ##   - in 'y': GL000225.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270729.1, KI270742.1, KI270754.1, KI270787.1, KI270803.1, ML143352.1, ML143375.1
    ##   - in 'y': GL000194.1, GL000216.2, GL000224.1, GL000254.2, KI270709.1, KI270732.1, KI270734.1, KI270738.1, KI270743.1, KI270770.1, KI270821.1, KI270924.1, KN196487.1, KN538361.1, KN538364.1, KN538366.1, KQ090027.1, KQ759759.1, KV575244.1, ML143341.1, ML143366.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270729.1, KI270742.1, KI270754.1, KI270787.1, KI270803.1, ML143352.1, ML143375.1
    ##   - in 'y': GL000194.1, GL000216.2, GL000224.1, GL000254.2, KI270709.1, KI270732.1, KI270734.1, KI270738.1, KI270743.1, KI270770.1, KI270821.1, KI270924.1, KN196487.1, KN538361.1, KN538364.1, KN538366.1, KQ090027.1, KQ759759.1, KV575244.1, ML143341.1, ML143366.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL383533.1, KI270442.1, KI270729.1, KI270742.1, KI270754.1, ML143352.1, KI270732.1, KI270738.1, KI270770.1, KN196487.1, KN538366.1, KQ759759.1, ML143372.1
    ##   - in 'y': GL000008.2, GL000225.1, KI270330.1, KI270717.1, KI270805.1, KI270880.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL383533.1, KI270442.1, KI270729.1, KI270742.1, KI270754.1, ML143352.1, KI270732.1, KI270738.1, KI270770.1, KN196487.1, KN538366.1, KQ759759.1, ML143372.1
    ##   - in 'y': GL000008.2, GL000225.1, KI270330.1, KI270717.1, KI270805.1, KI270880.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000251.2, KI270442.1, KI270742.1, KI270784.1, KI270787.1, KI270850.1, KI270902.1, ML143352.1, GL000216.2, KI270709.1, KI270732.1, KI270734.1, KI270738.1, KI270821.1, KI270924.1, KN196487.1, KN538361.1, KN538366.1, KQ090027.1, KQ759759.1, ML143341.1, ML143367.1, ML143372.1, GL000225.1, KI270330.1, KI270717.1, KI270805.1
    ##   - in 'y': KI270718.1, KI270849.1, KI270904.1, KN538370.1, KQ031389.1, KV880764.1, KZ559112.1, ML143353.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000251.2, KI270442.1, KI270742.1, KI270784.1, KI270787.1, KI270850.1, KI270902.1, ML143352.1, GL000216.2, KI270709.1, KI270732.1, KI270734.1, KI270738.1, KI270821.1, KI270924.1, KN196487.1, KN538361.1, KN538366.1, KQ090027.1, KQ759759.1, ML143341.1, ML143367.1, ML143372.1, GL000225.1, KI270330.1, KI270717.1, KI270805.1
    ##   - in 'y': KI270718.1, KI270849.1, KI270904.1, KN538370.1, KQ031389.1, KV880764.1, KZ559112.1, ML143353.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL383578.2, GL383581.2, JH159147.1, KI270712.1, KI270728.1, KI270765.1, KI270782.1, KI270783.1, KI270791.1, KI270803.1, KI270827.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270871.1, KI270878.1, KI270924.1, KN538364.1, KV766198.1, KV880764.1, KZ559105.1, ML143347.1, ML143353.1, ML143355.1, ML143359.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, GL383526.1, KI270311.1, KI270320.1, KI270330.1, KI270435.1, KI270589.1, KI270719.1, KI270724.1, KI270730.1, KI270732.1, KI270733.1, KI270736.1, KI270738.1, KI270754.1, KI270757.1, KI270787.1, KN196487.1, KQ031389.1, KQ759759.1, KQ983257.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL383578.2, GL383581.2, JH159147.1, KI270712.1, KI270728.1, KI270765.1, KI270782.1, KI270783.1, KI270791.1, KI270803.1, KI270827.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270871.1, KI270878.1, KI270924.1, KN538364.1, KV766198.1, KV880764.1, KZ559105.1, ML143347.1, ML143353.1, ML143355.1, ML143359.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, GL383526.1, KI270311.1, KI270320.1, KI270330.1, KI270435.1, KI270589.1, KI270719.1, KI270724.1, KI270730.1, KI270732.1, KI270733.1, KI270736.1, KI270738.1, KI270754.1, KI270757.1, KI270787.1, KN196487.1, KQ031389.1, KQ759759.1, KQ983257.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383581.2, KI270709.1, KI270779.1, KI270831.1, KI270836.1, KI270899.1, KN538370.1, KV880764.1, KZ208912.1, ML143350.1, ML143353.1, ML143365.1
    ##   - in 'y': GL000008.2, GL339449.2, KI270467.1, KI270707.1, KI270712.1, KI270714.1, KI270744.1, KI270830.1, KI270832.1, KI270849.1, KI270897.1, KI270908.1, KQ090028.1, KV766198.1, KZ559105.1, ML143367.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383581.2, KI270709.1, KI270779.1, KI270831.1, KI270836.1, KI270899.1, KN538370.1, KV880764.1, KZ208912.1, ML143350.1, ML143353.1, ML143365.1
    ##   - in 'y': GL000008.2, GL339449.2, KI270467.1, KI270707.1, KI270712.1, KI270714.1, KI270744.1, KI270830.1, KI270832.1, KI270849.1, KI270897.1, KI270908.1, KQ090028.1, KV766198.1, KZ559105.1, ML143367.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949746.1, KI270745.1, KI270856.1, KI270904.1, KI270908.1, KV575244.1, KV880768.1, KZ559105.1, ML143345.1
    ##   - in 'y': GL000008.2, KI270467.1, KI270728.1, KI270821.1, KQ031389.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949746.1, KI270745.1, KI270856.1, KI270904.1, KI270908.1, KV575244.1, KV880768.1, KZ559105.1, ML143345.1
    ##   - in 'y': GL000008.2, KI270467.1, KI270728.1, KI270821.1, KQ031389.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, JH159146.1, KI270337.1, KI270783.1, KI270836.1, KQ458384.1, KZ208907.1, chrY
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000254.2, GL383579.2, KI270438.1, KI270442.1, KI270711.1, KI270712.1, KI270721.1, KI270726.1, KI270733.1, KI270743.1, KI270745.1, KI270751.1, KI270754.1, KI270762.1, KI270804.1, KI270809.1, KI270816.1, KI270821.1, KI270851.1, KI270860.1, KI270861.1, KI270868.1, KI270880.1, KI270904.1, KN196487.1, KV880768.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383519.1, JH159146.1, KI270337.1, KI270783.1, KI270836.1, KQ458384.1, KZ208907.1, chrY
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000254.2, GL383579.2, KI270438.1, KI270442.1, KI270711.1, KI270712.1, KI270721.1, KI270726.1, KI270733.1, KI270743.1, KI270745.1, KI270751.1, KI270754.1, KI270762.1, KI270804.1, KI270809.1, KI270816.1, KI270821.1, KI270851.1, KI270860.1, KI270861.1, KI270868.1, KI270880.1, KI270904.1, KN196487.1, KV880768.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000253.2, GL383563.3, GL949752.1, JH159146.1, KI270438.1, KI270538.1, KI270709.1, KI270711.1, KI270728.1, KI270734.1, KI270760.1, KI270783.1, KI270784.1, KI270819.1, KI270831.1, KI270832.1, KI270850.1, KI270856.1, KI270868.1, KI270878.1, KI270903.1, KN196472.1, KQ458385.1, KV880763.1, KZ208912.1, KZ559112.1, ML143344.1, ML143366.1, ML143371.1, ML143377.1
    ##   - in 'y': KI270310.1, KI270706.1, KI270733.1, KI270849.1, KN196484.1, KV766198.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000253.2, GL383563.3, GL949752.1, JH159146.1, KI270438.1, KI270538.1, KI270709.1, KI270711.1, KI270728.1, KI270734.1, KI270760.1, KI270783.1, KI270784.1, KI270819.1, KI270831.1, KI270832.1, KI270850.1, KI270856.1, KI270868.1, KI270878.1, KI270903.1, KN196472.1, KQ458385.1, KV880763.1, KZ208912.1, KZ559112.1, ML143344.1, ML143366.1, ML143371.1, ML143377.1
    ##   - in 'y': KI270310.1, KI270706.1, KI270733.1, KI270849.1, KN196484.1, KV766198.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL949752.1, KI270714.1, KI270908.1, KV575244.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270730.1, KI270744.1, KI270757.1, KN196487.1, KQ983257.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000226.1, GL949752.1, KI270714.1, KI270908.1, KV575244.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270730.1, KI270744.1, KI270757.1, KN196487.1, KQ983257.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270742.1
    ##   - in 'y': GL000214.1, KI270707.1, KQ458385.1, chr13, chr16, chr18, chr20
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270742.1
    ##   - in 'y': GL000214.1, KI270707.1, KQ458385.1, chr13, chr16, chr18, chr20
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000251.2, GL383579.2, KI270438.1, KI270709.1, KI270714.1, KI270728.1, KI270765.1, KI270821.1, KI270848.1, KI270850.1, KI270861.1, KI270866.1, KI270897.1, KI270905.1, KV880764.1, KZ208913.1, KZ208921.1, ML143345.1
    ##   - in 'y': GL383578.2, GL949746.1, KI270733.1, KI270808.1, KI270849.1, KI270851.1, KI270892.1, KI270895.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000251.2, GL383579.2, KI270438.1, KI270709.1, KI270714.1, KI270728.1, KI270765.1, KI270821.1, KI270848.1, KI270850.1, KI270861.1, KI270866.1, KI270897.1, KI270905.1, KV880764.1, KZ208913.1, KZ208921.1, ML143345.1
    ##   - in 'y': GL383578.2, GL949746.1, KI270733.1, KI270808.1, KI270849.1, KI270851.1, KI270892.1, KI270895.1, KQ458383.1, KV880768.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000254.2, GL383580.2, KI270438.1, KI270821.1, KI270857.1, KI270861.1, KI270866.1, KI270897.1, KQ090026.1, KV766198.1, KZ208913.1, KZ208921.1, KI270733.1, KI270808.1, KI270849.1, KI270892.1, KI270895.1, KV880768.1, ML143375.1, chrM
    ##   - in 'y': GL000214.1, KI270738.1, KI270853.1, KI270856.1, KN538370.1, KQ031389.1, ML143352.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000254.2, GL383580.2, KI270438.1, KI270821.1, KI270857.1, KI270861.1, KI270866.1, KI270897.1, KQ090026.1, KV766198.1, KZ208913.1, KZ208921.1, KI270733.1, KI270808.1, KI270849.1, KI270892.1, KI270895.1, KV880768.1, ML143375.1, chrM
    ##   - in 'y': GL000214.1, KI270738.1, KI270853.1, KI270856.1, KN538370.1, KQ031389.1, ML143352.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000221.1, GL000251.2, GL000254.2, GL383579.2, GL383580.2, KI270438.1, KI270709.1, KI270714.1, KI270728.1, KI270821.1, KI270848.1, KI270850.1, KI270857.1, KI270861.1, KI270866.1, KI270897.1, KI270905.1, KN538364.1, KQ090026.1, KZ208913.1, KZ208921.1, ML143371.1, ML143380.1, GL383578.2, GL949746.1, KI270733.1, KI270808.1, KI270849.1, KI270851.1, KI270892.1, KI270895.1, KQ458383.1, KV880768.1, ML143375.1, chrM, GL000214.1, KI270738.1, KI270853.1, KI270856.1, ML143352.1, ML143365.1
    ##   - in 'y': KI270711.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000221.1, GL000251.2, GL000254.2, GL383579.2, GL383580.2, KI270438.1, KI270709.1, KI270714.1, KI270728.1, KI270821.1, KI270848.1, KI270850.1, KI270857.1, KI270861.1, KI270866.1, KI270897.1, KI270905.1, KN538364.1, KQ090026.1, KZ208913.1, KZ208921.1, ML143371.1, ML143380.1, GL383578.2, GL949746.1, KI270733.1, KI270808.1, KI270849.1, KI270851.1, KI270892.1, KI270895.1, KQ458383.1, KV880768.1, ML143375.1, chrM, GL000214.1, KI270738.1, KI270853.1, KI270856.1, ML143352.1, ML143365.1
    ##   - in 'y': KI270711.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL949746.1, GL949752.1, KI270721.1, KV575244.1, ML143345.1
    ##   - in 'y': GL000009.2, GL000218.1, GL383563.3, GL383578.2, KI270706.1, KI270783.1, KI270803.1, KI270844.1, KI270849.1, KI270908.1, KQ090026.1, KQ458383.1, KQ458384.1, KZ559112.1, ML143371.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL949746.1, GL949752.1, KI270721.1, KV575244.1, ML143345.1
    ##   - in 'y': GL000009.2, GL000218.1, GL383563.3, GL383578.2, KI270706.1, KI270783.1, KI270803.1, KI270844.1, KI270849.1, KI270908.1, KQ090026.1, KQ458383.1, KQ458384.1, KZ559112.1, ML143371.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270750.1, KI270860.1, KI270868.1, ML143345.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000214.1, GL000224.1, GL383581.2, KI270721.1, KI270722.1, KI270745.1, KI270767.1, KI270783.1, KI270803.1, KI270849.1, KI270856.1, KI270857.1, KI270899.1, KI270904.1, KN196484.1, KV575243.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270750.1, KI270860.1, KI270868.1, ML143345.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000214.1, GL000224.1, GL383581.2, KI270721.1, KI270722.1, KI270745.1, KI270767.1, KI270783.1, KI270803.1, KI270849.1, KI270856.1, KI270857.1, KI270899.1, KI270904.1, KN196484.1, KV575243.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000255.2, GL383578.2, KI270712.1, KI270728.1, KI270765.1, KI270849.1, KI270878.1, KN196479.1, KV766198.1, ML143344.1, ML143352.1, ML143355.1, ML143359.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, KI270310.1, KI270330.1, KI270435.1, KI270538.1, KI270590.1, KI270707.1, KI270709.1, KI270714.1, KI270729.1, KI270732.1, KI270736.1, KI270745.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270764.1, KI270783.1, KI270816.1, KI270819.1, KI270844.1, KI270851.1, KI270856.1, KI270908.1, KN196487.1, KQ031384.1, KQ983257.1, KV575244.1, KZ559112.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000255.2, GL383578.2, KI270712.1, KI270728.1, KI270765.1, KI270849.1, KI270878.1, KN196479.1, KV766198.1, ML143344.1, ML143352.1, ML143355.1, ML143359.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, KI270310.1, KI270330.1, KI270435.1, KI270538.1, KI270590.1, KI270707.1, KI270709.1, KI270714.1, KI270729.1, KI270732.1, KI270736.1, KI270745.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270764.1, KI270783.1, KI270816.1, KI270819.1, KI270844.1, KI270851.1, KI270856.1, KI270908.1, KN196487.1, KQ031384.1, KQ983257.1, KV575244.1, KZ559112.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KV766198.1, ML143380.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000221.1, GL000251.2, GL000255.2, GL949752.1, KI270438.1, KI270706.1, KI270711.1, KI270712.1, KI270714.1, KI270720.1, KI270721.1, KI270731.1, KI270744.1, KI270754.1, KI270763.1, KI270816.1, KI270819.1, KI270830.1, KI270844.1, KI270850.1, KI270853.1, KI270868.1, KI270879.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KQ458383.1, KV575244.1, KV880768.1, KZ559103.1, ML143371.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KV766198.1, ML143380.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000221.1, GL000251.2, GL000255.2, GL949752.1, KI270438.1, KI270706.1, KI270711.1, KI270712.1, KI270714.1, KI270720.1, KI270721.1, KI270731.1, KI270744.1, KI270754.1, KI270763.1, KI270816.1, KI270819.1, KI270830.1, KI270844.1, KI270850.1, KI270853.1, KI270868.1, KI270879.1, KI270905.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KQ458383.1, KV575244.1, KV880768.1, KZ559103.1, ML143371.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, GL383579.2, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270765.1, KI270783.1, KI270787.1, KI270832.1, KI270849.1, KI270861.1, KI270878.1, KI270899.1, KN538370.1, KQ031384.1, KQ458383.1, KV880764.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000221.1, GL000224.1, GL000225.1, GL339449.2, KI270709.1, KI270729.1, KI270731.1, KI270745.1, KI270764.1, KI270782.1, KI270831.1, KI270868.1, KI270869.1, KI270880.1, KI270897.1, KI270902.1, KI270903.1, KI270904.1, KN196484.1, KN196487.1, KN538372.1, KV766198.1, KZ208914.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, GL383579.2, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270765.1, KI270783.1, KI270787.1, KI270832.1, KI270849.1, KI270861.1, KI270878.1, KI270899.1, KN538370.1, KQ031384.1, KQ458383.1, KV880764.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000221.1, GL000224.1, GL000225.1, GL339449.2, KI270709.1, KI270729.1, KI270731.1, KI270745.1, KI270764.1, KI270782.1, KI270831.1, KI270868.1, KI270869.1, KI270880.1, KI270897.1, KI270902.1, KI270903.1, KI270904.1, KN196484.1, KN196487.1, KN538372.1, KV766198.1, KZ208914.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, KI270712.1, KI270718.1, KI270728.1, KI270765.1, KI270784.1, KI270857.1, KI270908.1, KN196479.1, KN538370.1, KV880764.1, ML143344.1, ML143352.1, ML143353.1, ML143359.1, ML143365.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000224.1, KI270442.1, KI270707.1, KI270711.1, KI270751.1, KI270782.1, KI270816.1, KI270844.1, KI270850.1, KI270851.1, KI270853.1, KI270897.1, KN196484.1, KQ090026.1, KV575244.1, KV880768.1, KZ559103.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383578.2, KI270712.1, KI270718.1, KI270728.1, KI270765.1, KI270784.1, KI270857.1, KI270908.1, KN196479.1, KN538370.1, KV880764.1, ML143344.1, ML143352.1, ML143353.1, ML143359.1, ML143365.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000224.1, KI270442.1, KI270707.1, KI270711.1, KI270751.1, KI270782.1, KI270816.1, KI270844.1, KI270850.1, KI270851.1, KI270853.1, KI270897.1, KN196484.1, KQ090026.1, KV575244.1, KV880768.1, KZ559103.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000253.2, GL339449.2, GL383563.3, GL383578.2, GL949746.1, GL949752.1, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270709.1, KI270711.1, KI270714.1, KI270726.1, KI270728.1, KI270731.1, KI270744.1, KI270745.1, KI270754.1, KI270770.1, KI270783.1, KI270784.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270819.1, KI270832.1, KI270844.1, KI270848.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270868.1, KI270878.1, KI270879.1, KN538364.1, KQ090026.1, KQ458383.1, KV766192.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143345.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270712.1, KI270729.1, KI270730.1, KI270897.1, ML143353.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000251.2, GL000253.2, GL339449.2, GL383563.3, GL383578.2, GL949746.1, GL949752.1, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270709.1, KI270711.1, KI270714.1, KI270726.1, KI270728.1, KI270731.1, KI270744.1, KI270745.1, KI270754.1, KI270770.1, KI270783.1, KI270784.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270819.1, KI270832.1, KI270844.1, KI270848.1, KI270849.1, KI270850.1, KI270853.1, KI270856.1, KI270857.1, KI270860.1, KI270861.1, KI270868.1, KI270878.1, KI270879.1, KN538364.1, KQ090026.1, KQ458383.1, KV766192.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143345.1, ML143366.1, ML143372.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270712.1, KI270729.1, KI270730.1, KI270897.1, ML143353.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270711.1, KI270731.1, KI270745.1, KI270787.1, KI270803.1, KI270832.1, KI270844.1, KI270878.1, KV766192.1, KZ208912.1, ML143366.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000221.1, GL000224.1, GL383522.1, GL383581.2, GL877875.1, KI270435.1, KI270438.1, KI270589.1, KI270718.1, KI270736.1, KI270738.1, KI270751.1, KI270753.1, KI270764.1, KI270765.1, KI270779.1, KI270782.1, KI270827.1, KI270831.1, KI270869.1, KI270880.1, KI270899.1, KI270902.1, KI270904.1, KI270938.1, KN196479.1, KN538360.1, KN538361.1, KN538370.1, KQ031388.1, KQ031389.1, KV880764.1, KZ208914.1, KZ208921.1, ML143344.1, ML143347.1, ML143350.1, ML143352.1, ML143355.1, ML143358.1, ML143362.1, ML143365.1, ML143370.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270711.1, KI270731.1, KI270745.1, KI270787.1, KI270803.1, KI270832.1, KI270844.1, KI270878.1, KV766192.1, KZ208912.1, ML143366.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000221.1, GL000224.1, GL383522.1, GL383581.2, GL877875.1, KI270435.1, KI270438.1, KI270589.1, KI270718.1, KI270736.1, KI270738.1, KI270751.1, KI270753.1, KI270764.1, KI270765.1, KI270779.1, KI270782.1, KI270827.1, KI270831.1, KI270869.1, KI270880.1, KI270899.1, KI270902.1, KI270904.1, KI270938.1, KN196479.1, KN538360.1, KN538361.1, KN538370.1, KQ031388.1, KQ031389.1, KV880764.1, KZ208914.1, KZ208921.1, ML143344.1, ML143347.1, ML143350.1, ML143352.1, ML143355.1, ML143358.1, ML143362.1, ML143365.1, ML143370.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000218.1, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270745.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270832.1, KI270844.1, KI270878.1, KV766192.1, KV766198.1, KZ208912.1, KZ559109.1, KI270897.1, GL000221.1, GL000224.1, GL383522.1, GL383581.2, GL877875.1, KI270435.1, KI270753.1, KI270827.1, KI270902.1, KI270904.1, KI270938.1, KN538361.1, KQ031388.1, KZ208914.1, ML143344.1, ML143347.1, ML143358.1, ML143362.1, ML143375.1
    ##   - in 'y': KI270330.1, KI270519.1, KI270591.1, KI270721.1, KI270851.1, KI270866.1, KI270937.1, KN538362.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000218.1, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270745.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270832.1, KI270844.1, KI270878.1, KV766192.1, KV766198.1, KZ208912.1, KZ559109.1, KI270897.1, GL000221.1, GL000224.1, GL383522.1, GL383581.2, GL877875.1, KI270435.1, KI270753.1, KI270827.1, KI270902.1, KI270904.1, KI270938.1, KN538361.1, KQ031388.1, KZ208914.1, ML143344.1, ML143347.1, ML143358.1, ML143362.1, ML143375.1
    ##   - in 'y': KI270330.1, KI270519.1, KI270591.1, KI270721.1, KI270851.1, KI270866.1, KI270937.1, KN538362.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270709.1, KI270711.1, KI270728.1, KI270733.1, KI270745.1, KI270754.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270844.1, KI270848.1, KV766192.1, KV766198.1, KZ208912.1, KI270712.1, KI270729.1, KI270730.1, ML143379.1, GL000008.2, GL000216.2, GL000221.1, GL000224.1, GL383522.1, GL383581.2, KI270435.1, KI270438.1, KI270589.1, KI270718.1, KI270736.1, KI270738.1, KI270751.1, KI270753.1, KI270764.1, KI270765.1, KI270779.1, KI270782.1, KI270827.1, KI270869.1, KI270880.1, KI270899.1, KI270902.1, KI270904.1, KI270938.1, KN196479.1, KN538360.1, KN538361.1, KN538370.1, KQ031388.1, KQ031389.1, KZ208914.1, KZ208921.1, ML143344.1, ML143347.1, ML143350.1, ML143352.1, ML143355.1, ML143358.1, ML143362.1, ML143365.1, ML143370.1, ML143375.1, KI270330.1, KI270519.1, KI270591.1, KI270721.1, KI270851.1, KI270866.1, KI270937.1, KN538362.1
    ##   - in 'y': GL000254.2, GL383520.2, JH159146.1, KQ458385.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL383578.2, KI270337.1, KI270466.1, KI270467.1, KI270707.1, KI270709.1, KI270711.1, KI270728.1, KI270733.1, KI270745.1, KI270754.1, KI270787.1, KI270803.1, KI270809.1, KI270816.1, KI270844.1, KI270848.1, KV766192.1, KV766198.1, KZ208912.1, KI270712.1, KI270729.1, KI270730.1, ML143379.1, GL000008.2, GL000216.2, GL000221.1, GL000224.1, GL383522.1, GL383581.2, KI270435.1, KI270438.1, KI270589.1, KI270718.1, KI270736.1, KI270738.1, KI270751.1, KI270753.1, KI270764.1, KI270765.1, KI270779.1, KI270782.1, KI270827.1, KI270869.1, KI270880.1, KI270899.1, KI270902.1, KI270904.1, KI270938.1, KN196479.1, KN538360.1, KN538361.1, KN538370.1, KQ031388.1, KQ031389.1, KZ208914.1, KZ208921.1, ML143344.1, ML143347.1, ML143350.1, ML143352.1, ML143355.1, ML143358.1, ML143362.1, ML143365.1, ML143370.1, ML143375.1, KI270330.1, KI270519.1, KI270591.1, KI270721.1, KI270851.1, KI270866.1, KI270937.1, KN538362.1
    ##   - in 'y': GL000254.2, GL383520.2, JH159146.1, KQ458385.1, KZ208922.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000214.1, GL000216.2, GL000220.1, GL383578.2, GL383581.2, KI270330.1, KI270435.1, KI270709.1, KI270712.1, KI270721.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270751.1, KI270764.1, KI270765.1, KI270782.1, KI270869.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KN538370.1, KQ031384.1, KQ031389.1, ML143344.1, ML143345.1, ML143347.1, ML143352.1, ML143353.1, ML143355.1, ML143358.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000218.1, GL339449.2, GL383563.3, KI270707.1, KI270714.1, KI270745.1, KI270754.1, KI270762.1, KI270787.1, KI270816.1, KI270850.1, KI270856.1, KI270861.1, KI270876.1, KI270878.1, KI270894.1, KI270897.1, KN196484.1, KV880768.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000214.1, GL000216.2, GL000220.1, GL383578.2, GL383581.2, KI270330.1, KI270435.1, KI270709.1, KI270712.1, KI270721.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270751.1, KI270764.1, KI270765.1, KI270782.1, KI270869.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KN538370.1, KQ031384.1, KQ031389.1, ML143344.1, ML143345.1, ML143347.1, ML143352.1, ML143353.1, ML143355.1, ML143358.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000218.1, GL339449.2, GL383563.3, KI270707.1, KI270714.1, KI270745.1, KI270754.1, KI270762.1, KI270787.1, KI270816.1, KI270850.1, KI270856.1, KI270861.1, KI270876.1, KI270878.1, KI270894.1, KI270897.1, KN196484.1, KV880768.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270435.1, KI270712.1, KI270738.1, KI270751.1, KI270764.1, KI270765.1, KI270831.1, KI270936.1, KI270937.1, KQ031384.1, ML143359.1, KI270762.1, KI270876.1, KI270878.1, KN196484.1
    ##   - in 'y': GL000251.2, GL000254.2, GL383522.1, GL877875.1, GL949752.1, KI270538.1, KI270711.1, KI270717.1, KI270726.1, KI270730.1, KI270731.1, KI270770.1, KI270779.1, KI270784.1, KI270791.1, KI270792.1, KI270803.1, KI270804.1, KI270807.1, KI270813.1, KI270821.1, KI270827.1, KI270832.1, KI270835.1, KI270844.1, KI270848.1, KI270857.1, KI270868.1, KI270872.1, KI270892.1, KI270895.1, KI270903.1, KI270938.1, KN196479.1, KN196480.1, KN538364.1, KV766198.1, KV880763.1, KZ208908.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, KZ559105.1, ML143350.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270435.1, KI270712.1, KI270738.1, KI270751.1, KI270764.1, KI270765.1, KI270831.1, KI270936.1, KI270937.1, KQ031384.1, ML143359.1, KI270762.1, KI270876.1, KI270878.1, KN196484.1
    ##   - in 'y': GL000251.2, GL000254.2, GL383522.1, GL877875.1, GL949752.1, KI270538.1, KI270711.1, KI270717.1, KI270726.1, KI270730.1, KI270731.1, KI270770.1, KI270779.1, KI270784.1, KI270791.1, KI270792.1, KI270803.1, KI270804.1, KI270807.1, KI270813.1, KI270821.1, KI270827.1, KI270832.1, KI270835.1, KI270844.1, KI270848.1, KI270857.1, KI270868.1, KI270872.1, KI270892.1, KI270895.1, KI270903.1, KI270938.1, KN196479.1, KN196480.1, KN538364.1, KV766198.1, KV880763.1, KZ208908.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, KZ559105.1, ML143350.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000216.2, GL000220.1, GL000224.1, KI270330.1, KI270435.1, KI270438.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270751.1, KI270764.1, KI270869.1, KI270899.1, KI270936.1, KI270937.1, KN538370.1, KQ031384.1, ML143344.1, ML143347.1, ML143353.1, ML143358.1, ML143359.1, ML143366.1, ML143375.1, ML143379.1, GL000218.1, KI270787.1, KI270816.1, KI270876.1, KI270894.1, KN196484.1, GL000254.2, GL877875.1, KI270538.1, KI270717.1, KI270726.1, KI270730.1, KI270731.1, KI270791.1, KI270792.1, KI270803.1, KI270804.1, KI270813.1, KI270827.1, KI270832.1, KI270835.1, KI270844.1, KI270872.1, KI270892.1, KI270895.1, KI270903.1, KI270938.1, KN196479.1, KV880763.1, KZ208908.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, KZ559105.1, ML143350.1, ML143365.1, ML143370.1
    ##   - in 'y': GL383580.2, JH159146.1, KI270734.1, KI270809.1, KI270904.1, KN538372.1, KQ090016.1, KV766193.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000194.1, GL000216.2, GL000220.1, GL000224.1, KI270330.1, KI270435.1, KI270438.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270751.1, KI270764.1, KI270869.1, KI270899.1, KI270936.1, KI270937.1, KN538370.1, KQ031384.1, ML143344.1, ML143347.1, ML143353.1, ML143358.1, ML143359.1, ML143366.1, ML143375.1, ML143379.1, GL000218.1, KI270787.1, KI270816.1, KI270876.1, KI270894.1, KN196484.1, GL000254.2, GL877875.1, KI270538.1, KI270717.1, KI270726.1, KI270730.1, KI270731.1, KI270791.1, KI270792.1, KI270803.1, KI270804.1, KI270813.1, KI270827.1, KI270832.1, KI270835.1, KI270844.1, KI270872.1, KI270892.1, KI270895.1, KI270903.1, KI270938.1, KN196479.1, KV880763.1, KZ208908.1, KZ208912.1, KZ208914.1, KZ208921.1, KZ559103.1, KZ559105.1, ML143350.1, ML143365.1, ML143370.1
    ##   - in 'y': GL383580.2, JH159146.1, KI270734.1, KI270809.1, KI270904.1, KN538372.1, KQ090016.1, KV766193.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KQ031384.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000221.1, GL000225.1, GL000253.2, GL000255.2, GL383563.3, KI270538.1, KI270714.1, KI270742.1, KI270754.1, KI270783.1, KI270792.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270868.1, KI270905.1, KI270908.1, KN538364.1, KV575244.1, KV766198.1, KZ208912.1, KZ208915.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270438.1, KQ031384.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000221.1, GL000225.1, GL000253.2, GL000255.2, GL383563.3, KI270538.1, KI270714.1, KI270742.1, KI270754.1, KI270783.1, KI270792.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270868.1, KI270905.1, KI270908.1, KN538364.1, KV575244.1, KV766198.1, KZ208912.1, KZ208915.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, KI270734.1, KI270765.1, KI270816.1, KI270832.1, KI270844.1, KI270848.1, KI270868.1, KI270937.1, KN538361.1, KQ458384.1, KV880764.1, KZ559105.1, KZ559109.1, ML143377.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL000251.2, GL383578.2, KI270707.1, KI270712.1, KI270728.1, KI270729.1, KI270730.1, KI270849.1, KI270861.1, KI270934.1, KN538369.1, KQ031389.1, KQ458383.1, KZ559112.1, ML143352.1, ML143353.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, KI270734.1, KI270765.1, KI270816.1, KI270832.1, KI270844.1, KI270848.1, KI270868.1, KI270937.1, KN538361.1, KQ458384.1, KV880764.1, KZ559105.1, KZ559109.1, ML143377.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL000251.2, GL383578.2, KI270707.1, KI270712.1, KI270728.1, KI270729.1, KI270730.1, KI270849.1, KI270861.1, KI270934.1, KN538369.1, KQ031389.1, KQ458383.1, KZ559112.1, ML143352.1, ML143353.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000253.2, KI270734.1, KI270744.1, KI270848.1, KZ559105.1, KZ559109.1, GL000194.1, GL000224.1, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270849.1, KI270934.1, KQ458383.1, chrM
    ##   - in 'y': GL383574.1, KI270718.1, KI270721.1, KI270731.1, KI270806.1, KI270865.1, KI270908.1, KV880768.1, KZ208921.1, ML143359.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000253.2, KI270734.1, KI270744.1, KI270848.1, KZ559105.1, KZ559109.1, GL000194.1, GL000224.1, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270849.1, KI270934.1, KQ458383.1, chrM
    ##   - in 'y': GL383574.1, KI270718.1, KI270721.1, KI270731.1, KI270806.1, KI270865.1, KI270908.1, KV880768.1, KZ208921.1, ML143359.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL339449.2, GL383522.1, GL383581.2, GL877875.1, KI270734.1, KI270745.1, KI270784.1, KI270816.1, KI270830.1, KI270832.1, KI270844.1, KI270848.1, KI270866.1, KN196479.1, KN538361.1, KQ458384.1, KV575244.1, KV766198.1, KZ559105.1, KZ559109.1, ML143366.1, ML143380.1, GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL383578.2, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270849.1, KI270934.1, KN538369.1, KZ559112.1, ML143375.1, chrM, KI270721.1, KI270731.1, KI270806.1, KI270865.1, KV880768.1, KZ208921.1
    ##   - in 'y': GL000216.2, KI270438.1, KI270709.1, KI270720.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL339449.2, GL383522.1, GL383581.2, GL877875.1, KI270734.1, KI270745.1, KI270784.1, KI270816.1, KI270830.1, KI270832.1, KI270844.1, KI270848.1, KI270866.1, KN196479.1, KN538361.1, KQ458384.1, KV575244.1, KV766198.1, KZ559105.1, KZ559109.1, ML143366.1, ML143380.1, GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL383578.2, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270849.1, KI270934.1, KN538369.1, KZ559112.1, ML143375.1, chrM, KI270721.1, KI270731.1, KI270806.1, KI270865.1, KV880768.1, KZ208921.1
    ##   - in 'y': GL000216.2, KI270438.1, KI270709.1, KI270720.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270830.1, KI270853.1, KQ458385.1, ML143345.1
    ##   - in 'y': GL000255.2, GL339449.2, KI270728.1, KI270849.1, KQ458383.1, KV766198.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, KI270830.1, KI270853.1, KQ458385.1, ML143345.1
    ##   - in 'y': GL000255.2, GL339449.2, KI270728.1, KI270849.1, KQ458383.1, KV766198.1, ML143372.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000251.2, GL000256.2, GL383578.2, KI270337.1, KI270466.1, KI270712.1, KI270718.1, KI270728.1, KI270733.1, KI270735.1, KI270742.1, KI270750.1, KI270751.1, KI270765.1, KI270803.1, KI270857.1, KI270892.1, KI270894.1, KI270899.1, KN196479.1, KN196487.1, KN538364.1, KQ031384.1, KV575244.1, KV766198.1, KV880764.1, KZ208915.1, ML143345.1, ML143352.1, ML143355.1, ML143359.1, ML143365.1, ML143366.1, ML143371.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000009.2, GL000216.2, GL000220.1, GL000224.1, KI270442.1, KI270754.1, KI270850.1, KI270851.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000251.2, GL000256.2, GL383578.2, KI270337.1, KI270466.1, KI270712.1, KI270718.1, KI270728.1, KI270733.1, KI270735.1, KI270742.1, KI270750.1, KI270751.1, KI270765.1, KI270803.1, KI270857.1, KI270892.1, KI270894.1, KI270899.1, KN196479.1, KN196487.1, KN538364.1, KQ031384.1, KV575244.1, KV766198.1, KV880764.1, KZ208915.1, ML143345.1, ML143352.1, ML143355.1, ML143359.1, ML143365.1, ML143366.1, ML143371.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000009.2, GL000216.2, GL000220.1, GL000224.1, KI270442.1, KI270754.1, KI270850.1, KI270851.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000257.2, GL383556.1, GL877875.1, KI270720.1, KI270730.1, KI270751.1, KI270802.1, KI270810.1, KI270838.1, KI270862.1, KI270938.1, KN196472.1, KV766192.1, KV766193.1, KV880763.1, KV880765.1, KZ559103.1, ML143343.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000221.1, GL339449.2, GL949752.1, KI270712.1, KI270719.1, KI270731.1, KI270733.1, KI270745.1, KI270765.1, KI270782.1, KI270804.1, KI270819.1, KI270821.1, KI270831.1, KI270849.1, KI270866.1, KI270880.1, KI270899.1, KI270904.1, KI270937.1, KQ031389.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143367.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000257.2, GL383556.1, GL877875.1, KI270720.1, KI270730.1, KI270751.1, KI270802.1, KI270810.1, KI270838.1, KI270862.1, KI270938.1, KN196472.1, KV766192.1, KV766193.1, KV880763.1, KV880765.1, KZ559103.1, ML143343.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000221.1, GL339449.2, GL949752.1, KI270712.1, KI270719.1, KI270731.1, KI270733.1, KI270745.1, KI270765.1, KI270782.1, KI270804.1, KI270819.1, KI270821.1, KI270831.1, KI270849.1, KI270866.1, KI270880.1, KI270899.1, KI270904.1, KI270937.1, KQ031389.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143367.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000257.2, GL383522.1, GL383556.1, KI270707.1, KI270720.1, KI270799.1, KI270802.1, KI270810.1, KI270838.1, KN196472.1, KN538372.1, KQ458385.1, KV766192.1, KV880763.1, KV880765.1, KZ208913.1, KZ559103.1, ML143343.1, GL000194.1, GL000218.1, GL000221.1, GL949752.1, KI270712.1, KI270731.1, KI270733.1, KI270831.1, KI270849.1, KI270866.1, KI270880.1, KI270899.1, KI270937.1, ML143344.1, ML143350.1, ML143365.1, ML143375.1
    ##   - in 'y': GL383526.1, GL383557.1, GL383581.2, KI270591.1, KI270709.1, KI270734.1, KI270735.1, KI270738.1, KI270754.1, KI270779.1, KI270806.1, KI270854.1, KI270861.1, KI270924.1, KN538360.1, KQ458386.1, KV575243.1, ML143354.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000257.2, GL383522.1, GL383556.1, KI270707.1, KI270720.1, KI270799.1, KI270802.1, KI270810.1, KI270838.1, KN196472.1, KN538372.1, KQ458385.1, KV766192.1, KV880763.1, KV880765.1, KZ208913.1, KZ559103.1, ML143343.1, GL000194.1, GL000218.1, GL000221.1, GL949752.1, KI270712.1, KI270731.1, KI270733.1, KI270831.1, KI270849.1, KI270866.1, KI270880.1, KI270899.1, KI270937.1, ML143344.1, ML143350.1, ML143365.1, ML143375.1
    ##   - in 'y': GL383526.1, GL383557.1, GL383581.2, KI270591.1, KI270709.1, KI270734.1, KI270735.1, KI270738.1, KI270754.1, KI270779.1, KI270806.1, KI270854.1, KI270861.1, KI270924.1, KN538360.1, KQ458386.1, KV575243.1, ML143354.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, KI270438.1, KI270832.1, KI270838.1, KI270862.1, KI270938.1, KQ458384.1, KV880763.1, KZ559103.1, GL000194.1, GL000218.1, GL000221.1, KI270712.1, KI270731.1, KI270733.1, KI270765.1, KI270819.1, KI270821.1, KI270831.1, KI270866.1, KI270880.1, KI270899.1, KI270937.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143367.1, ML143375.1, ML143379.1, GL383526.1, GL383557.1, GL383581.2, KI270591.1, KI270709.1, KI270734.1, KI270735.1, KI270738.1, KI270754.1, KI270779.1, KI270806.1, KI270854.1, KI270861.1, KI270924.1, KN538360.1, KQ458386.1, KV575243.1, ML143354.1, ML143370.1
    ##   - in 'y': GL000252.2, GL383546.1, KI270589.1, KI270784.1, KI270844.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, KI270438.1, KI270832.1, KI270838.1, KI270862.1, KI270938.1, KQ458384.1, KV880763.1, KZ559103.1, GL000194.1, GL000218.1, GL000221.1, KI270712.1, KI270731.1, KI270733.1, KI270765.1, KI270819.1, KI270821.1, KI270831.1, KI270866.1, KI270880.1, KI270899.1, KI270937.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143367.1, ML143375.1, ML143379.1, GL383526.1, GL383557.1, GL383581.2, KI270591.1, KI270709.1, KI270734.1, KI270735.1, KI270738.1, KI270754.1, KI270779.1, KI270806.1, KI270854.1, KI270861.1, KI270924.1, KN538360.1, KQ458386.1, KV575243.1, ML143354.1, ML143370.1
    ##   - in 'y': GL000252.2, GL383546.1, KI270589.1, KI270784.1, KI270844.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, KI270733.1, KI270744.1, ML143366.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000224.1, GL000225.1, KI270709.1, KI270742.1, KI270908.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, KI270733.1, KI270744.1, ML143366.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000224.1, GL000225.1, KI270709.1, KI270742.1, KI270908.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000258.2, GL383567.1, KI270330.1, KI270582.1, KI270741.1, KI270745.1, KI270810.1, KI270871.1, KN196480.1, KN538370.1, KQ458383.1, KZ208912.1, ML143352.1, ML143353.1, ML143354.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000254.2, KI270707.1, KI270709.1, KI270712.1, KI270717.1, KI270734.1, KI270754.1, KI270762.1, KI270763.1, KI270779.1, KI270781.1, KI270783.1, KI270787.1, KI270816.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270869.1, KI270870.1, KI270903.1, KI270926.1, KI270937.1, KN196484.1, KN196487.1, KQ983257.1, KV880763.1, KZ559103.1, ML143360.1, ML143364.1, ML143366.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000258.2, GL383567.1, KI270330.1, KI270582.1, KI270741.1, KI270745.1, KI270810.1, KI270871.1, KN196480.1, KN538370.1, KQ458383.1, KZ208912.1, ML143352.1, ML143353.1, ML143354.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000254.2, KI270707.1, KI270709.1, KI270712.1, KI270717.1, KI270734.1, KI270754.1, KI270762.1, KI270763.1, KI270779.1, KI270781.1, KI270783.1, KI270787.1, KI270816.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270869.1, KI270870.1, KI270903.1, KI270926.1, KI270937.1, KN196484.1, KN196487.1, KQ983257.1, KV880763.1, KZ559103.1, ML143360.1, ML143364.1, ML143366.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000258.2, GL383526.1, GL383567.1, GL877875.1, KI270330.1, KI270442.1, KI270582.1, KI270718.1, KI270720.1, KI270721.1, KI270728.1, KI270741.1, KI270745.1, KI270784.1, KI270803.1, KI270810.1, KI270830.1, KI270868.1, KI270871.1, KI270880.1, KN196480.1, KQ458383.1, KQ458384.1, KZ208912.1, KZ208921.1, ML143354.1, GL000254.2, KI270734.1, KI270762.1, KI270763.1, KI270781.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270870.1, KI270903.1, KI270926.1, KI270937.1, KN196484.1, KN196487.1, KQ983257.1, KZ559103.1, ML143360.1, ML143364.1, ML143366.1, ML143380.1
    ##   - in 'y': GL383581.2, KI270310.1, KI270731.1, KI270736.1, KI270821.1, KI270848.1, KI270892.1, KI270924.1, KI270938.1, KN538364.1, KQ458385.1, KV880768.1, ML143347.1, ML143358.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000258.2, GL383526.1, GL383567.1, GL877875.1, KI270330.1, KI270442.1, KI270582.1, KI270718.1, KI270720.1, KI270721.1, KI270728.1, KI270741.1, KI270745.1, KI270784.1, KI270803.1, KI270810.1, KI270830.1, KI270868.1, KI270871.1, KI270880.1, KN196480.1, KQ458383.1, KQ458384.1, KZ208912.1, KZ208921.1, ML143354.1, GL000254.2, KI270734.1, KI270762.1, KI270763.1, KI270781.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270870.1, KI270903.1, KI270926.1, KI270937.1, KN196484.1, KN196487.1, KQ983257.1, KZ559103.1, ML143360.1, ML143364.1, ML143366.1, ML143380.1
    ##   - in 'y': GL383581.2, KI270310.1, KI270731.1, KI270736.1, KI270821.1, KI270848.1, KI270892.1, KI270924.1, KI270938.1, KN538364.1, KQ458385.1, KV880768.1, ML143347.1, ML143358.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000218.1, GL000220.1, GL000251.2, GL000258.2, GL383522.1, GL383567.1, GL383578.2, KI270330.1, KI270438.1, KI270582.1, KI270720.1, KI270733.1, KI270741.1, KI270745.1, KI270767.1, KI270803.1, KI270804.1, KI270810.1, KI270868.1, KI270880.1, KN196480.1, KQ458383.1, KQ458384.1, KZ208912.1, KZ208921.1, ML143354.1, ML143372.1, GL000225.1, GL000254.2, KI270707.1, KI270709.1, KI270762.1, KI270763.1, KI270781.1, KI270783.1, KI270816.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270870.1, KI270903.1, KI270926.1, KN196484.1, KN196487.1, KQ983257.1, KZ559103.1, ML143364.1, ML143366.1, ML143380.1, GL383581.2, KI270310.1, KI270731.1, KI270736.1, KI270821.1, KI270848.1, KI270892.1, KI270924.1, KQ458385.1, ML143347.1
    ##   - in 'y': KI270711.1, KI270726.1, KI270730.1, KI270899.1, KQ031384.1, KQ031389.1, KZ208913.1, ML143344.1, ML143350.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000218.1, GL000220.1, GL000251.2, GL000258.2, GL383522.1, GL383567.1, GL383578.2, KI270330.1, KI270438.1, KI270582.1, KI270720.1, KI270733.1, KI270741.1, KI270745.1, KI270767.1, KI270803.1, KI270804.1, KI270810.1, KI270868.1, KI270880.1, KN196480.1, KQ458383.1, KQ458384.1, KZ208912.1, KZ208921.1, ML143354.1, ML143372.1, GL000225.1, GL000254.2, KI270707.1, KI270709.1, KI270762.1, KI270763.1, KI270781.1, KI270783.1, KI270816.1, KI270818.1, KI270844.1, KI270851.1, KI270860.1, KI270865.1, KI270870.1, KI270903.1, KI270926.1, KN196484.1, KN196487.1, KQ983257.1, KZ559103.1, ML143364.1, ML143366.1, ML143380.1, GL383581.2, KI270310.1, KI270731.1, KI270736.1, KI270821.1, KI270848.1, KI270892.1, KI270924.1, KQ458385.1, ML143347.1
    ##   - in 'y': KI270711.1, KI270726.1, KI270730.1, KI270899.1, KQ031384.1, KQ031389.1, KZ208913.1, ML143344.1, ML143350.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, KI270337.1, KI270742.1, KI270765.1, KV575244.1, ML143352.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000225.1, KI270438.1, KI270442.1, KI270709.1, KI270744.1, KI270751.1, KN196487.1, KQ090026.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000256.2, KI270337.1, KI270742.1, KI270765.1, KV575244.1, ML143352.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000225.1, KI270438.1, KI270442.1, KI270709.1, KI270744.1, KI270751.1, KN196487.1, KQ090026.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000251.2, GL000254.2, GL383519.1, GL383542.1, GL383567.1, GL949752.1, KI270721.1, KI270726.1, KI270743.1, KI270770.1, KI270784.1, KI270792.1, KI270827.1, KI270844.1, KI270860.1, KI270863.1, KI270879.1, KI270903.1, KN196484.1, KN538361.1, KQ031389.1, KV766192.1, KV880764.1, KZ208921.1, KZ559105.1, KZ559109.1, ML143345.1, ML143353.1
    ##   - in 'y': GL000220.1, GL383527.1, GL383578.2, GL383580.2, KI270310.1, KI270322.1, KI270842.1, KI270896.1, KZ208913.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000205.2, GL000251.2, GL000254.2, GL383519.1, GL383542.1, GL383567.1, GL949752.1, KI270721.1, KI270726.1, KI270743.1, KI270770.1, KI270784.1, KI270792.1, KI270827.1, KI270844.1, KI270860.1, KI270863.1, KI270879.1, KI270903.1, KN196484.1, KN538361.1, KQ031389.1, KV766192.1, KV880764.1, KZ208921.1, KZ559105.1, KZ559109.1, ML143345.1, ML143353.1
    ##   - in 'y': GL000220.1, GL383527.1, GL383578.2, GL383580.2, KI270310.1, KI270322.1, KI270842.1, KI270896.1, KZ208913.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270857.1, KV880764.1, ML143345.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, GL000225.1, GL383563.3, GL877875.1, KI270709.1, KI270736.1, KI270742.1, KI270744.1, KI270751.1, KI270754.1, KI270782.1, KI270836.1, KI270880.1, KN196484.1, KN538364.1, KQ031389.1, ML143365.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270857.1, KV880764.1, ML143345.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, GL000225.1, GL383563.3, GL877875.1, KI270709.1, KI270736.1, KI270742.1, KI270744.1, KI270751.1, KI270754.1, KI270782.1, KI270836.1, KI270880.1, KN196484.1, KN538364.1, KQ031389.1, ML143365.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL383563.3, KI270442.1, KI270706.1, KI270707.1, KI270709.1, KI270714.1, KI270728.1, KI270742.1, KI270765.1, KI270783.1, KI270784.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270861.1, KI270905.1, KI270908.1, KQ090026.1, KV766198.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, ML143345.1, ML143353.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270337.1, KI270467.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL383563.3, KI270442.1, KI270706.1, KI270707.1, KI270709.1, KI270714.1, KI270728.1, KI270742.1, KI270765.1, KI270783.1, KI270784.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270861.1, KI270905.1, KI270908.1, KQ090026.1, KV766198.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, ML143345.1, ML143353.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   - in 'y': KI270337.1, KI270467.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL383563.3, KI270442.1, KI270706.1, KI270707.1, KI270709.1, KI270714.1, KI270744.1, KI270783.1, KI270784.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270861.1, KI270905.1, KI270908.1, KQ090026.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, ML143353.1, ML143372.1, ML143377.1, ML143380.1, chrM, KI270337.1, KI270467.1
    ##   - in 'y': KI270712.1, KZ559109.1, ML143355.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, GL383563.3, KI270442.1, KI270706.1, KI270707.1, KI270709.1, KI270714.1, KI270744.1, KI270783.1, KI270784.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270861.1, KI270905.1, KI270908.1, KQ090026.1, KV880768.1, KZ208912.1, KZ208915.1, KZ559105.1, ML143353.1, ML143372.1, ML143377.1, ML143380.1, chrM, KI270337.1, KI270467.1
    ##   - in 'y': KI270712.1, KZ559109.1, ML143355.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270707.1, KI270709.1, KQ090026.1, KV880768.1, KZ208912.1, KZ559105.1, ML143353.1, KI270337.1, KI270467.1, KZ559109.1, ML143375.1
    ##   - in 'y': GL000205.2, GL000251.2, KI270745.1, KI270832.1, KI270879.1, KN538364.1, KN538370.1, KQ031389.1, KV880764.1, KZ559112.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270707.1, KI270709.1, KQ090026.1, KV880768.1, KZ208912.1, KZ559105.1, ML143353.1, KI270337.1, KI270467.1, KZ559109.1, ML143375.1
    ##   - in 'y': GL000205.2, GL000251.2, KI270745.1, KI270832.1, KI270879.1, KN538364.1, KN538370.1, KQ031389.1, KV880764.1, KZ559112.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270934.1, KI270936.1, KI270937.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000225.1, KI270320.1, KI270330.1, KI270435.1, KI270442.1, KI270465.1, KI270508.1, KI270510.1, KI270512.1, KI270515.1, KI270519.1, KI270538.1, KI270589.1, KI270590.1, KI270591.1, KI270712.1, KI270719.1, KI270723.1, KI270728.1, KI270729.1, KI270730.1, KI270732.1, KI270733.1, KI270735.1, KI270736.1, KI270737.1, KI270738.1, KI270742.1, KI270743.1, KI270744.1, KI270746.1, KI270751.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270772.1, KI270894.1, KN196487.1, KN538372.1, KQ031384.1, KQ983257.1, ML143354.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1, chr10, chr13, chr14, chr18, chr20, chr21, chr6
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270934.1, KI270936.1, KI270937.1
    ##   - in 'y': GL000008.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000225.1, KI270320.1, KI270330.1, KI270435.1, KI270442.1, KI270465.1, KI270508.1, KI270510.1, KI270512.1, KI270515.1, KI270519.1, KI270538.1, KI270589.1, KI270590.1, KI270591.1, KI270712.1, KI270719.1, KI270723.1, KI270728.1, KI270729.1, KI270730.1, KI270732.1, KI270733.1, KI270735.1, KI270736.1, KI270737.1, KI270738.1, KI270742.1, KI270743.1, KI270744.1, KI270746.1, KI270751.1, KI270754.1, KI270756.1, KI270757.1, KI270762.1, KI270772.1, KI270894.1, KN196487.1, KN538372.1, KQ031384.1, KQ983257.1, ML143354.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1, chr10, chr13, chr14, chr18, chr20, chr21, chr6
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270934.1, KI270936.1, KI270937.1, KI270465.1, KI270510.1, KI270515.1, KI270590.1, KI270737.1, ML143371.1, ML143372.1, ML143377.1
    ##   - in 'y': GL000009.2, GL000208.1, KI270310.1, KI270713.1, KI270714.1, KI270718.1, KI270722.1, KI270724.1, KI270725.1, KI270731.1, KI270750.1, KI270770.1, KI270779.1, KI270805.1, KI270836.1, KI270857.1, KI270880.1, KN538364.1, KQ031389.1, KQ759759.1, KZ208913.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270934.1, KI270936.1, KI270937.1, KI270465.1, KI270510.1, KI270515.1, KI270590.1, KI270737.1, ML143371.1, ML143372.1, ML143377.1
    ##   - in 'y': GL000009.2, GL000208.1, KI270310.1, KI270713.1, KI270714.1, KI270718.1, KI270722.1, KI270724.1, KI270725.1, KI270731.1, KI270750.1, KI270770.1, KI270779.1, KI270805.1, KI270836.1, KI270857.1, KI270880.1, KN538364.1, KQ031389.1, KQ759759.1, KZ208913.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270718.1, KI270770.1, KI270861.1, KI270905.1, KI270936.1, KN538360.1, KQ759759.1, KV766192.1, ML143344.1, ML143355.1
    ##   - in 'y': GL000194.1, KI270310.1, KI270311.1, KI270538.1, KI270712.1, KI270713.1, KI270728.1, KI270731.1, KI270746.1, KI270757.1, KI270762.1, KI270805.1, KI270836.1, KI270880.1, KI270894.1, KN196472.1, KN196487.1, KN538372.1, KQ458385.1, KQ983257.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270718.1, KI270770.1, KI270861.1, KI270905.1, KI270936.1, KN538360.1, KQ759759.1, KV766192.1, ML143344.1, ML143355.1
    ##   - in 'y': GL000194.1, KI270310.1, KI270311.1, KI270538.1, KI270712.1, KI270713.1, KI270728.1, KI270731.1, KI270746.1, KI270757.1, KI270762.1, KI270805.1, KI270836.1, KI270880.1, KI270894.1, KN196472.1, KN196487.1, KN538372.1, KQ458385.1, KQ983257.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, KI270465.1, KI270508.1, KI270519.1, KI270538.1, KI270587.1, KI270712.1, KI270713.1, KI270714.1, KI270724.1, KI270731.1, KI270732.1, KI270738.1, KI270894.1, KN538372.1, KQ983257.1, KZ208915.1, ML143380.1, chrX
    ##   - in 'y': KI270718.1, KI270762.1, KI270782.1, KQ031384.1, KV880764.1, ML143354.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000194.1, KI270465.1, KI270508.1, KI270519.1, KI270538.1, KI270587.1, KI270712.1, KI270713.1, KI270714.1, KI270724.1, KI270731.1, KI270732.1, KI270738.1, KI270894.1, KN538372.1, KQ983257.1, KZ208915.1, ML143380.1, chrX
    ##   - in 'y': KI270718.1, KI270762.1, KI270782.1, KQ031384.1, KV880764.1, ML143354.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': ML143352.1, ML143354.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270713.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, chr14
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': ML143352.1, ML143354.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270713.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, chr14
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KZ208915.1, ML143352.1, ML143354.1, KI270465.1, KI270853.1, KI270908.1, KN196487.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270591.1, KI270712.1, KI270714.1, KI270723.1, KI270724.1, KI270728.1, KI270729.1, KI270730.1, KI270731.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270742.1, KI270756.1, KI270762.1, KI270805.1, KI270836.1, KI270894.1, KI270902.1, KN538364.1, KN538372.1, KQ031384.1, KQ983257.1, KZ208913.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KZ208915.1, ML143352.1, ML143354.1, KI270465.1, KI270853.1, KI270908.1, KN196487.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270591.1, KI270712.1, KI270714.1, KI270723.1, KI270724.1, KI270728.1, KI270729.1, KI270730.1, KI270731.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270742.1, KI270756.1, KI270762.1, KI270805.1, KI270836.1, KI270894.1, KI270902.1, KN538364.1, KN538372.1, KQ031384.1, KQ983257.1, KZ208913.1, ML143372.1, ML143377.1, ML143378.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270733.1, KZ208915.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270723.1, KI270724.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270756.1, KI270805.1, KI270836.1, KI270894.1, KN538364.1, KN538372.1, KQ983257.1, KZ208913.1, ML143372.1, ML143377.1, ML143378.1
    ##   - in 'y': KI270718.1, KI270770.1, KI270779.1, KI270808.1, KI270895.1, KI270899.1, KV880764.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270733.1, KZ208915.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270723.1, KI270724.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270756.1, KI270805.1, KI270836.1, KI270894.1, KN538364.1, KN538372.1, KQ983257.1, KZ208913.1, ML143372.1, ML143377.1, ML143378.1
    ##   - in 'y': KI270718.1, KI270770.1, KI270779.1, KI270808.1, KI270895.1, KI270899.1, KV880764.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270330.1, KI270733.1, KI270736.1, KI270744.1, KZ208915.1, ML143354.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270713.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270591.1, KI270714.1, KI270723.1, KI270724.1, KI270729.1, KI270730.1, KI270731.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270742.1, KI270756.1, KI270762.1, KI270805.1, KI270836.1, KI270902.1, KN538364.1, KQ031384.1, KQ983257.1, KZ208913.1, ML143377.1, ML143378.1, KI270770.1, KI270779.1, KI270808.1, KI270895.1, KI270899.1, KV880764.1, ML143355.1, ML143379.1
    ##   - in 'y': KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270330.1, KI270733.1, KI270736.1, KI270744.1, KZ208915.1, ML143354.1, GL000225.1, KI270435.1, KI270442.1, KI270465.1, KI270713.1, KI270751.1, KI270757.1, KI270853.1, KI270908.1, KN196487.1, GL000009.2, GL000194.1, GL000219.1, KI270310.1, KI270508.1, KI270519.1, KI270538.1, KI270589.1, KI270591.1, KI270714.1, KI270723.1, KI270724.1, KI270729.1, KI270730.1, KI270731.1, KI270732.1, KI270734.1, KI270735.1, KI270738.1, KI270742.1, KI270756.1, KI270762.1, KI270805.1, KI270836.1, KI270902.1, KN538364.1, KQ031384.1, KQ983257.1, KZ208913.1, ML143377.1, ML143378.1, KI270770.1, KI270779.1, KI270808.1, KI270895.1, KI270899.1, KV880764.1, ML143355.1, ML143379.1
    ##   - in 'y': KI270754.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, KI270538.1, KI270827.1, KQ983257.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383526.1, GL383578.2, GL877875.1, KI270467.1, KI270719.1, KI270720.1, KI270731.1, KI270732.1, KI270743.1, KI270751.1, KI270803.1, KI270817.1, KI270821.1, KI270849.1, KI270850.1, KI270878.1, KI270895.1, KI270905.1, KN196487.1, KN538364.1, KQ458385.1, KQ983255.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, KI270538.1, KI270827.1, KQ983257.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383526.1, GL383578.2, GL877875.1, KI270467.1, KI270719.1, KI270720.1, KI270731.1, KI270732.1, KI270743.1, KI270751.1, KI270803.1, KI270817.1, KI270821.1, KI270849.1, KI270850.1, KI270878.1, KI270895.1, KI270905.1, KN196487.1, KN538364.1, KQ458385.1, KQ983255.1, KZ559111.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL339449.2, KI270330.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270538.1, KI270707.1, KI270709.1, KI270712.1, KI270718.1, KI270720.1, KI270729.1, KI270730.1, KI270735.1, KI270737.1, KI270738.1, KI270751.1, KI270765.1, KI270842.1, KI270909.1, KQ031384.1, KZ559100.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000225.1, GL000254.2, GL949752.1, KI270714.1, KI270721.1, KI270783.1, KI270792.1, KI270802.1, KI270830.1, KI270832.1, KI270844.1, KI270849.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270903.1, KI270904.1, KN196484.1, KN196487.1, KQ090026.1, KZ208913.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL339449.2, KI270330.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270538.1, KI270707.1, KI270709.1, KI270712.1, KI270718.1, KI270720.1, KI270729.1, KI270730.1, KI270735.1, KI270737.1, KI270738.1, KI270751.1, KI270765.1, KI270842.1, KI270909.1, KQ031384.1, KZ559100.1, ML143344.1, ML143350.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000225.1, GL000254.2, GL949752.1, KI270714.1, KI270721.1, KI270783.1, KI270792.1, KI270802.1, KI270830.1, KI270832.1, KI270844.1, KI270849.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270903.1, KI270904.1, KN196484.1, KN196487.1, KQ090026.1, KZ208913.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270734.1, KI270765.1, KI270842.1, KN538361.1, ML143355.1
    ##   - in 'y': GL000214.1, GL000224.1, GL339449.2, GL383581.2, GL877875.1, GL949746.1, GL949752.1, KI270310.1, KI270411.1, KI270442.1, KI270466.1, KI270467.1, KI270711.1, KI270712.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270784.1, KI270821.1, KI270830.1, KI270849.1, KI270853.1, KI270871.1, KI270892.1, KI270897.1, KI270908.1, KN538364.1, KV575244.1, KZ208912.1, KZ208915.1, ML143350.1, ML143353.1, ML143358.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143375.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270734.1, KI270765.1, KI270842.1, KN538361.1, ML143355.1
    ##   - in 'y': GL000214.1, GL000224.1, GL339449.2, GL383581.2, GL877875.1, GL949746.1, GL949752.1, KI270310.1, KI270411.1, KI270442.1, KI270466.1, KI270467.1, KI270711.1, KI270712.1, KI270744.1, KI270745.1, KI270754.1, KI270782.1, KI270784.1, KI270821.1, KI270830.1, KI270849.1, KI270853.1, KI270871.1, KI270892.1, KI270897.1, KI270908.1, KN538364.1, KV575244.1, KZ208912.1, KZ208915.1, ML143350.1, ML143353.1, ML143358.1, ML143365.1, ML143366.1, ML143371.1, ML143372.1, ML143375.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270466.1, KI270850.1, KN196476.1, ML143352.1, ML143373.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000252.2, GL383556.1, GL949752.1, KI270438.1, KI270711.1, KI270713.1, KI270714.1, KI270728.1, KI270738.1, KI270753.1, KI270787.1, KI270809.1, KI270822.1, KI270856.1, KI270868.1, KI270876.1, KI270901.1, KI270905.1, KI270924.1, KN538364.1, KQ458383.1, KQ458384.1, KZ208913.1, KZ559109.1, ML143344.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270466.1, KI270850.1, KN196476.1, ML143352.1, ML143373.1
    ##   - in 'y': GL000214.1, GL000225.1, GL000252.2, GL383556.1, GL949752.1, KI270438.1, KI270711.1, KI270713.1, KI270714.1, KI270728.1, KI270738.1, KI270753.1, KI270787.1, KI270809.1, KI270822.1, KI270856.1, KI270868.1, KI270876.1, KI270901.1, KI270905.1, KI270924.1, KN538364.1, KQ458383.1, KQ458384.1, KZ208913.1, KZ559109.1, ML143344.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270709.1, KI270751.1, KI270754.1, KI270784.1, KI270819.1, KI270866.1, KI270895.1
    ##   - in 'y': GL000221.1, GL000251.2, GL000254.2, GL383563.3, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270742.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270856.1, KI270865.1, KI270868.1, KI270869.1, KI270871.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270899.1, KI270937.1, KN196479.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143341.1, ML143358.1, ML143359.1, ML143362.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270709.1, KI270751.1, KI270754.1, KI270784.1, KI270819.1, KI270866.1, KI270895.1
    ##   - in 'y': GL000221.1, GL000251.2, GL000254.2, GL383563.3, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270742.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270856.1, KI270865.1, KI270868.1, KI270869.1, KI270871.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270899.1, KI270937.1, KN196479.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143341.1, ML143358.1, ML143359.1, ML143362.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270714.1, KI270754.1, KI270784.1, KI270819.1, KI270895.1, ML143344.1, ML143375.1, GL000221.1, GL000251.2, GL000254.2, GL383563.3, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270742.1, KI270783.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270851.1, KI270865.1, KI270868.1, KI270869.1, KI270871.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270937.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143341.1, ML143359.1, ML143362.1, ML143364.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   - in 'y': KI270330.1, KI270538.1, KI270772.1, KI270907.1, KQ031384.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270714.1, KI270754.1, KI270784.1, KI270819.1, KI270895.1, ML143344.1, ML143375.1, GL000221.1, GL000251.2, GL000254.2, GL383563.3, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270742.1, KI270783.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270849.1, KI270850.1, KI270851.1, KI270865.1, KI270868.1, KI270869.1, KI270871.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270937.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143341.1, ML143359.1, ML143362.1, ML143364.1, ML143372.1, ML143377.1, ML143380.1, chrM
    ##   - in 'y': KI270330.1, KI270538.1, KI270772.1, KI270907.1, KQ031384.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL339449.2, GL949746.1, KI270709.1, KI270718.1, KI270733.1, KI270736.1, KI270751.1, KI270784.1, KI270819.1, KI270895.1, KN538360.1, GL000221.1, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270830.1, KI270831.1, KI270832.1, KI270865.1, KI270869.1, KI270871.1, KI270880.1, KI270897.1, KI270899.1, KI270937.1, KN196479.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143358.1, ML143364.1, ML143372.1, ML143380.1, chrM, KI270330.1, KI270538.1, KI270772.1, KI270907.1, KQ031384.1
    ##   - in 'y': KI270734.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL339449.2, GL949746.1, KI270709.1, KI270718.1, KI270733.1, KI270736.1, KI270751.1, KI270784.1, KI270819.1, KI270895.1, KN538360.1, GL000221.1, GL383567.1, GL383578.2, KI270442.1, KI270589.1, KI270707.1, KI270711.1, KI270730.1, KI270830.1, KI270831.1, KI270832.1, KI270865.1, KI270869.1, KI270871.1, KI270880.1, KI270897.1, KI270899.1, KI270937.1, KN196479.1, KQ090026.1, KQ458383.1, KZ559109.1, ML143358.1, ML143364.1, ML143372.1, ML143380.1, chrM, KI270330.1, KI270538.1, KI270772.1, KI270907.1, KQ031384.1
    ##   - in 'y': KI270734.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000254.2, GL339449.2, GL383563.3, KI270337.1, KI270467.1, KI270714.1, KI270805.1, KI270809.1, KI270878.1, KI270880.1, KI270904.1, KI270907.1, KI270924.1, KQ031389.1, KZ559112.1, ML143352.1, ML143353.1, ML143364.1, ML143366.1, ML143379.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270706.1, KI270709.1, KI270711.1, KI270726.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KQ090026.1, KZ208912.1, KZ208914.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000254.2, GL339449.2, GL383563.3, KI270337.1, KI270467.1, KI270714.1, KI270805.1, KI270809.1, KI270878.1, KI270880.1, KI270904.1, KI270907.1, KI270924.1, KQ031389.1, KZ559112.1, ML143352.1, ML143353.1, ML143364.1, ML143366.1, ML143379.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270706.1, KI270709.1, KI270711.1, KI270726.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KQ090026.1, KZ208912.1, KZ208914.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, KI270337.1, KI270442.1, KI270467.1, KI270830.1, KI270848.1, KI270857.1, KI270860.1, KI270878.1, KI270904.1, KI270924.1, KQ031389.1, KV766198.1, KZ559112.1, ML143352.1, ML143353.1, ML143364.1, ML143366.1, ML143379.1, GL000216.2, GL000224.1, KI270519.1, KI270706.1, KI270711.1, KI270726.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KZ208912.1, KZ208914.1
    ##   - in 'y': KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, KI270337.1, KI270442.1, KI270467.1, KI270830.1, KI270848.1, KI270857.1, KI270860.1, KI270878.1, KI270904.1, KI270924.1, KQ031389.1, KV766198.1, KZ559112.1, ML143352.1, ML143353.1, ML143364.1, ML143366.1, ML143379.1, GL000216.2, GL000224.1, KI270519.1, KI270706.1, KI270711.1, KI270726.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KZ208912.1, KZ208914.1
    ##   - in 'y': KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270337.1, KI270467.1, KI270804.1, KI270805.1, KI270904.1, KI270907.1, KI270924.1, ML143352.1, ML143364.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270709.1, KI270711.1, KI270844.1, KN196487.1, KZ208912.1, KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1
    ##   - in 'y': GL383580.2, GL383581.2, GL877875.1, KI270707.1, KI270721.1, KI270728.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KI270803.1, KI270821.1, KI270832.1, KI270897.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270337.1, KI270467.1, KI270804.1, KI270805.1, KI270904.1, KI270907.1, KI270924.1, ML143352.1, ML143364.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270709.1, KI270711.1, KI270844.1, KN196487.1, KZ208912.1, KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1
    ##   - in 'y': GL383580.2, GL383581.2, GL877875.1, KI270707.1, KI270721.1, KI270728.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KI270803.1, KI270821.1, KI270832.1, KI270897.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270337.1, KI270467.1, KI270805.1, KI270848.1, KI270878.1, KI270904.1, KI270907.1, KI270924.1, KQ031389.1, ML143352.1, ML143353.1, ML143364.1, ML143371.1, ML143379.1, GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270709.1, KI270844.1, KN538364.1, KZ208912.1, KI270330.1, KI270724.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1, ML143345.1, GL383580.2, GL383581.2, GL877875.1, KI270707.1, KI270721.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1
    ##   - in 'y': GL383578.2, KI270827.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, KI270337.1, KI270467.1, KI270805.1, KI270848.1, KI270878.1, KI270904.1, KI270907.1, KI270924.1, KQ031389.1, ML143352.1, ML143353.1, ML143364.1, ML143371.1, ML143379.1, GL000214.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270709.1, KI270844.1, KN538364.1, KZ208912.1, KI270330.1, KI270724.1, KI270733.1, KI270761.1, KI270836.1, KI270892.1, ML143345.1, GL383580.2, GL383581.2, GL877875.1, KI270707.1, KI270721.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1
    ##   - in 'y': GL383578.2, KI270827.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000254.2, GL383563.3, KI270337.1, KI270442.1, KI270467.1, KI270804.1, KI270805.1, KI270809.1, KI270831.1, KI270878.1, KI270880.1, KI270904.1, KI270905.1, KI270907.1, KI270924.1, KV766198.1, ML143352.1, ML143364.1, ML143366.1, ML143367.1, ML143379.1, GL000194.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270706.1, KI270709.1, KI270711.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KQ090026.1, KZ208912.1, KZ208914.1, KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, GL383580.2, GL383581.2, KI270707.1, KI270721.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KI270803.1, KI270821.1, KI270832.1, KI270897.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1, GL383578.2, KI270827.1
    ##   - in 'y': KI270720.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000254.2, GL383563.3, KI270337.1, KI270442.1, KI270467.1, KI270804.1, KI270805.1, KI270809.1, KI270831.1, KI270878.1, KI270880.1, KI270904.1, KI270905.1, KI270907.1, KI270924.1, KV766198.1, ML143352.1, ML143364.1, ML143366.1, ML143367.1, ML143379.1, GL000194.1, GL000216.2, GL000220.1, GL000221.1, GL000224.1, GL000225.1, KI270438.1, KI270519.1, KI270706.1, KI270709.1, KI270711.1, KI270816.1, KI270844.1, KN196487.1, KN538364.1, KQ090026.1, KZ208912.1, KZ208914.1, KI270330.1, KI270724.1, KI270731.1, KI270733.1, KI270761.1, KI270836.1, GL383580.2, GL383581.2, KI270707.1, KI270721.1, KI270734.1, KI270762.1, KI270765.1, KI270792.1, KI270803.1, KI270821.1, KI270832.1, KI270897.1, KN538361.1, KQ458383.1, KV880764.1, KZ559105.1, ML143355.1, ML143375.1, GL383578.2, KI270827.1
    ##   - in 'y': KI270720.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270735.1, KZ208908.1, ML143352.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000216.2, GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL000255.2, GL383578.2, GL383583.2, GL949746.1, KI270442.1, KI270707.1, KI270745.1, KI270779.1, KI270782.1, KI270816.1, KI270819.1, KI270831.1, KI270844.1, KI270849.1, KI270853.1, KI270879.1, KI270880.1, KI270905.1, KI270908.1, KI270924.1, KI270937.1, KN196484.1, KN196487.1, KN538364.1, KQ090026.1, KQ458383.1, KV575244.1, KV880764.1, KZ208912.1, KZ208922.1, KZ559112.1, ML143367.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270735.1, KZ208908.1, ML143352.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000216.2, GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL000255.2, GL383578.2, GL383583.2, GL949746.1, KI270442.1, KI270707.1, KI270745.1, KI270779.1, KI270782.1, KI270816.1, KI270819.1, KI270831.1, KI270844.1, KI270849.1, KI270853.1, KI270879.1, KI270880.1, KI270905.1, KI270908.1, KI270924.1, KI270937.1, KN196484.1, KN196487.1, KN538364.1, KQ090026.1, KQ458383.1, KV575244.1, KV880764.1, KZ208912.1, KZ208922.1, KZ559112.1, ML143367.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000221.1, GL383522.1, GL949746.1, GL949752.1, KI270438.1, KI270709.1, KI270734.1, KI270736.1, KI270738.1, KI270754.1, KI270764.1, KI270779.1, KI270844.1, KI270848.1, KI270850.1, KI270856.1, KI270866.1, KI270868.1, KN538370.1, KQ031389.1, KV575244.1, KV880764.1, KV880768.1, KZ208914.1, KZ559105.1, ML143345.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143366.1, ML143375.1, ML143380.1
    ##   - in 'y': GL000009.2, KI270712.1, KI270783.1, KI270803.1, KI270813.1, KI270851.1, KI270865.1, KI270879.1, KI270899.1, KQ090026.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000221.1, GL383522.1, GL949746.1, GL949752.1, KI270438.1, KI270709.1, KI270734.1, KI270736.1, KI270738.1, KI270754.1, KI270764.1, KI270779.1, KI270844.1, KI270848.1, KI270850.1, KI270856.1, KI270866.1, KI270868.1, KN538370.1, KQ031389.1, KV575244.1, KV880764.1, KV880768.1, KZ208914.1, KZ559105.1, ML143345.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143366.1, ML143375.1, ML143380.1
    ##   - in 'y': GL000009.2, KI270712.1, KI270783.1, KI270803.1, KI270813.1, KI270851.1, KI270865.1, KI270879.1, KI270899.1, KQ090026.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, GL949752.1, KI270438.1, KI270442.1, KI270728.1, KI270733.1, KI270734.1, KI270736.1, KI270738.1, KI270754.1, KI270764.1, KI270765.1, KI270779.1, KI270807.1, KI270844.1, KI270848.1, KI270849.1, KI270850.1, KI270856.1, KI270857.1, KI270866.1, KN538364.1, KN538370.1, KQ458383.1, KV766198.1, KV880764.1, KV880768.1, KZ208914.1, KZ559105.1, ML143352.1, ML143353.1, ML143365.1, ML143366.1, ML143375.1, ML143380.1, KI270712.1, KI270851.1, KI270865.1, KI270879.1
    ##   - in 'y': GL000251.2, GL000254.2, KZ208922.1, ML143344.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, GL949752.1, KI270438.1, KI270442.1, KI270728.1, KI270733.1, KI270734.1, KI270736.1, KI270738.1, KI270754.1, KI270764.1, KI270765.1, KI270779.1, KI270807.1, KI270844.1, KI270848.1, KI270849.1, KI270850.1, KI270856.1, KI270857.1, KI270866.1, KN538364.1, KN538370.1, KQ458383.1, KV766198.1, KV880764.1, KV880768.1, KZ208914.1, KZ559105.1, ML143352.1, ML143353.1, ML143365.1, ML143366.1, ML143375.1, ML143380.1, KI270712.1, KI270851.1, KI270865.1, KI270879.1
    ##   - in 'y': GL000251.2, GL000254.2, KZ208922.1, ML143344.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL949752.1, KI270721.1, KI270736.1, KI270856.1, KI270866.1, KZ208914.1, KZ559105.1, ML143365.1, ML143366.1, GL000009.2, KI270712.1, KI270783.1, KI270813.1, KI270851.1, KI270865.1, KI270879.1, KQ090026.1, ML143372.1, GL000251.2, GL000254.2, KZ208922.1, ML143344.1
    ##   - in 'y': GL000008.2, GL000216.2, GL339449.2, KI270718.1, KI270730.1, KI270735.1, KI270745.1, KI270782.1, KI270784.1, KI270819.1, KI270831.1, KI270869.1, KI270924.1, KI270936.1, KN538360.1, KN538361.1, KN538362.1, ML143350.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL949752.1, KI270721.1, KI270736.1, KI270856.1, KI270866.1, KZ208914.1, KZ559105.1, ML143365.1, ML143366.1, GL000009.2, KI270712.1, KI270783.1, KI270813.1, KI270851.1, KI270865.1, KI270879.1, KQ090026.1, ML143372.1, GL000251.2, GL000254.2, KZ208922.1, ML143344.1
    ##   - in 'y': GL000008.2, GL000216.2, GL339449.2, KI270718.1, KI270730.1, KI270735.1, KI270745.1, KI270782.1, KI270784.1, KI270819.1, KI270831.1, KI270869.1, KI270924.1, KI270936.1, KN538360.1, KN538361.1, KN538362.1, ML143350.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000205.2, GL000214.1, GL000218.1, GL000224.1, GL000225.1, GL383578.2, KI270310.1, KI270438.1, KI270538.1, KI270707.1, KI270709.1, KI270712.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270743.1, KI270751.1, KI270754.1, KI270762.1, KI270770.1, KI270779.1, KI270816.1, KI270832.1, KI270861.1, KI270866.1, KI270868.1, KI270878.1, KI270892.1, KI270903.1, KN196480.1, KN538361.1, KN538372.1, KQ031384.1, KQ983257.1, KV766193.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143366.1, ML143377.1
    ##   - in 'y': GL000009.2, GL339449.2, GL383581.2, KI270337.1, KI270467.1, KI270821.1, KI270844.1, KV880763.1, KZ559103.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000205.2, GL000214.1, GL000218.1, GL000224.1, GL000225.1, GL383578.2, KI270310.1, KI270438.1, KI270538.1, KI270707.1, KI270709.1, KI270712.1, KI270719.1, KI270721.1, KI270728.1, KI270733.1, KI270734.1, KI270743.1, KI270751.1, KI270754.1, KI270762.1, KI270770.1, KI270779.1, KI270816.1, KI270832.1, KI270861.1, KI270866.1, KI270868.1, KI270878.1, KI270892.1, KI270903.1, KN196480.1, KN538361.1, KN538372.1, KQ031384.1, KQ983257.1, KV766193.1, KV880768.1, KZ208915.1, KZ208922.1, KZ559105.1, ML143366.1, ML143377.1
    ##   - in 'y': GL000009.2, GL339449.2, GL383581.2, KI270337.1, KI270467.1, KI270821.1, KI270844.1, KV880763.1, KZ559103.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000251.2, KI270337.1, KI270467.1, KI270849.1, KI270879.1, KI270908.1, KV766198.1, KZ208915.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270330.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270744.1, KI270754.1, KI270850.1, KN196487.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000251.2, KI270337.1, KI270467.1, KI270849.1, KI270879.1, KI270908.1, KV766198.1, KZ208915.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270330.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270736.1, KI270738.1, KI270744.1, KI270754.1, KI270850.1, KN196487.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000254.2, GL383522.1, GL949746.1, KI270761.1, KI270787.1, KI270850.1, KI270868.1, KQ090026.1, KV575244.1, KV880764.1, ML143369.1
    ##   - in 'y': GL339449.2, GL383526.1, GL383580.2, KI270742.1, KI270783.1, KI270784.1, KI270821.1, KI270892.1, KQ031389.1, KQ458383.1, ML143366.1, ML143372.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000254.2, GL383522.1, GL949746.1, KI270761.1, KI270787.1, KI270850.1, KI270868.1, KQ090026.1, KV575244.1, KV880764.1, ML143369.1
    ##   - in 'y': GL339449.2, GL383526.1, GL383580.2, KI270742.1, KI270783.1, KI270784.1, KI270821.1, KI270892.1, KQ031389.1, KQ458383.1, ML143366.1, ML143372.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949752.1, KI270442.1, KI270721.1, KI270728.1, KI270742.1, KI270782.1, KI270783.1, KI270830.1, KI270861.1, KI270879.1, KI270897.1, KN538364.1, KN538370.1, KV880764.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000219.1, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143371.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949752.1, KI270442.1, KI270721.1, KI270728.1, KI270742.1, KI270782.1, KI270783.1, KI270830.1, KI270861.1, KI270879.1, KI270897.1, KN538364.1, KN538370.1, KV880764.1, ML143353.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000219.1, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143371.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383578.2, GL949752.1, KI270711.1, KI270721.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270861.1, KI270897.1, KN538364.1, KN538370.1, KQ458383.1, KV880764.1, ML143353.1, ML143375.1, chrM, GL000219.1, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143372.1
    ##   - in 'y': KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, KV880768.1, ML143344.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383578.2, GL949752.1, KI270711.1, KI270721.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270861.1, KI270897.1, KN538364.1, KN538370.1, KQ458383.1, KV880764.1, ML143353.1, ML143375.1, chrM, GL000219.1, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143372.1
    ##   - in 'y': KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, KV880768.1, ML143344.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270744.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270853.1, KI270861.1, KI270868.1, KI270879.1, KI270897.1, KI270908.1, KN538364.1, KQ090026.1, KQ458383.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143345.1, ML143353.1, ML143355.1, ML143375.1, ML143380.1, chrM, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143371.1, ML143372.1, KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, KV880768.1, ML143344.1, ML143377.1
    ##   - in 'y': KN538362.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270721.1, KI270744.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270853.1, KI270861.1, KI270868.1, KI270879.1, KI270897.1, KI270908.1, KN538364.1, KQ090026.1, KQ458383.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143345.1, ML143353.1, ML143355.1, ML143375.1, ML143380.1, chrM, GL949746.1, KI270337.1, KI270466.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, KV575244.1, ML143371.1, ML143372.1, KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, KV880768.1, ML143344.1, ML143377.1
    ##   - in 'y': KN538362.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383578.2, GL949752.1, KI270442.1, KI270706.1, KI270711.1, KI270721.1, KI270742.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270853.1, KI270861.1, KI270868.1, KI270897.1, KI270908.1, KN538364.1, KQ458383.1, KV766198.1, KZ559112.1, ML143353.1, ML143355.1, ML143375.1, chrM, GL000219.1, GL949746.1, KI270337.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, ML143371.1, ML143372.1, KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, ML143344.1, ML143377.1, KN538362.1
    ##   - in 'y': GL000214.1, KI270438.1, KI270718.1, KI270733.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383578.2, GL949752.1, KI270442.1, KI270706.1, KI270711.1, KI270721.1, KI270742.1, KI270782.1, KI270783.1, KI270784.1, KI270830.1, KI270853.1, KI270861.1, KI270868.1, KI270897.1, KI270908.1, KN538364.1, KQ458383.1, KV766198.1, KZ559112.1, ML143353.1, ML143355.1, ML143375.1, chrM, GL000219.1, GL949746.1, KI270337.1, KI270467.1, KI270731.1, KI270792.1, KI270813.1, KI270849.1, KI270856.1, KI270857.1, KI270860.1, KI270880.1, KI270899.1, KI270936.1, KI270937.1, KQ031389.1, KQ458384.1, ML143371.1, ML143372.1, KI270745.1, KI270765.1, KI270770.1, KI270779.1, KI270819.1, ML143344.1, ML143377.1, KN538362.1
    ##   - in 'y': GL000214.1, KI270438.1, KI270718.1, KI270733.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270333.1, KI270442.1, KI270711.1, KI270850.1, KI270894.1, KN196487.1, KQ031384.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000219.1, GL000220.1, GL000224.1, KI270742.1, KI270744.1, KI270754.1, KQ031389.1, KQ983257.1, ML143371.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270333.1, KI270442.1, KI270711.1, KI270850.1, KI270894.1, KN196487.1, KQ031384.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000219.1, GL000220.1, GL000224.1, KI270742.1, KI270744.1, KI270754.1, KQ031389.1, KQ983257.1, ML143371.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270333.1, KI270711.1, KI270713.1, GL000009.2, GL000194.1, KI270742.1, KQ031389.1, KQ983257.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000225.1, KI270330.1, KI270435.1, KI270465.1, KI270508.1, KI270538.1, KI270591.1, KI270706.1, KI270707.1, KI270712.1, KI270728.1, KI270729.1, KI270730.1, KI270731.1, KI270733.1, KI270735.1, KI270736.1, KI270751.1, KI270757.1, KI270762.1, KI270899.1, KN538364.1, KN538372.1, ML143345.1, ML143353.1, ML143356.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270333.1, KI270711.1, KI270713.1, GL000009.2, GL000194.1, KI270742.1, KQ031389.1, KQ983257.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000225.1, KI270330.1, KI270435.1, KI270465.1, KI270508.1, KI270538.1, KI270591.1, KI270706.1, KI270707.1, KI270712.1, KI270728.1, KI270729.1, KI270730.1, KI270731.1, KI270733.1, KI270735.1, KI270736.1, KI270751.1, KI270757.1, KI270762.1, KI270899.1, KN538364.1, KN538372.1, ML143345.1, ML143353.1, ML143356.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, KI270333.1, KI270442.1, KI270466.1, KI270467.1, KI270709.1, KI270711.1, KI270713.1, KI270850.1, KI270894.1, KN196487.1, KQ031384.1, GL000009.2, GL000194.1, GL000219.1, GL000220.1, KI270744.1, KQ983257.1, ML143371.1, ML143377.1, GL000008.2, GL000216.2, GL000225.1, KI270435.1, KI270465.1, KI270508.1, KI270538.1, KI270591.1, KI270706.1, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270731.1, KI270733.1, KI270735.1, KI270736.1, KI270751.1, KI270757.1, KI270762.1, KI270899.1, KN538364.1, KN538372.1, ML143356.1, ML143372.1
    ##   - in 'y': KI270718.1, KI270908.1, KQ090026.1, KV880764.1, KZ559112.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, KI270333.1, KI270442.1, KI270466.1, KI270467.1, KI270709.1, KI270711.1, KI270713.1, KI270850.1, KI270894.1, KN196487.1, KQ031384.1, GL000009.2, GL000194.1, GL000219.1, GL000220.1, KI270744.1, KQ983257.1, ML143371.1, ML143377.1, GL000008.2, GL000216.2, GL000225.1, KI270435.1, KI270465.1, KI270508.1, KI270538.1, KI270591.1, KI270706.1, KI270707.1, KI270712.1, KI270729.1, KI270730.1, KI270731.1, KI270733.1, KI270735.1, KI270736.1, KI270751.1, KI270757.1, KI270762.1, KI270899.1, KN538364.1, KN538372.1, ML143356.1, ML143372.1
    ##   - in 'y': KI270718.1, KI270908.1, KQ090026.1, KV880764.1, KZ559112.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270310.1, KI270337.1, KN538372.1, KQ031384.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000225.1, GL000255.2, KI270442.1, KI270706.1, KI270712.1, KI270736.1, KI270743.1, KI270744.1, KI270757.1, KI270816.1, KI270850.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270310.1, KI270337.1, KN538372.1, KQ031384.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000225.1, GL000255.2, KI270442.1, KI270706.1, KI270712.1, KI270736.1, KI270743.1, KI270744.1, KI270757.1, KI270816.1, KI270850.1, KZ559105.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, ML143352.1, ML143359.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000253.2, GL000255.2, GL949752.1, KI270330.1, KI270442.1, KI270707.1, KI270709.1, KI270714.1, KI270721.1, KI270726.1, KI270729.1, KI270733.1, KI270742.1, KI270751.1, KI270754.1, KI270764.1, KI270783.1, KI270813.1, KI270816.1, KI270819.1, KI270830.1, KI270832.1, KI270844.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270857.1, KI270878.1, KI270879.1, KI270880.1, KI270894.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KV575244.1, KV766198.1, KZ208915.1, KZ559112.1, ML143353.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, ML143382.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, ML143352.1, ML143359.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000251.2, GL000253.2, GL000255.2, GL949752.1, KI270330.1, KI270442.1, KI270707.1, KI270709.1, KI270714.1, KI270721.1, KI270726.1, KI270729.1, KI270733.1, KI270742.1, KI270751.1, KI270754.1, KI270764.1, KI270783.1, KI270813.1, KI270816.1, KI270819.1, KI270830.1, KI270832.1, KI270844.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270857.1, KI270878.1, KI270879.1, KI270880.1, KI270894.1, KI270908.1, KN196484.1, KN196487.1, KN538364.1, KV575244.1, KV766198.1, KZ208915.1, KZ559112.1, ML143353.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1, ML143382.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270712.1, ML143352.1, ML143359.1, GL000194.1, GL000220.1, GL000224.1, GL000225.1, GL949752.1, KI270330.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270754.1, KI270764.1, KI270816.1, KI270819.1, KI270832.1, KI270894.1, KN196484.1, KN196487.1, KN538364.1, ML143366.1
    ##   - in 'y': GL000216.2, KI270587.1, KI270711.1, KI270720.1, KI270745.1, KI270803.1, KI270892.1, KI270897.1, KQ458383.1, KQ458385.1, ML143355.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270712.1, ML143352.1, ML143359.1, GL000194.1, GL000220.1, GL000224.1, GL000225.1, GL949752.1, KI270330.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270754.1, KI270764.1, KI270816.1, KI270819.1, KI270832.1, KI270894.1, KN196484.1, KN196487.1, KN538364.1, ML143366.1
    ##   - in 'y': GL000216.2, KI270587.1, KI270711.1, KI270720.1, KI270745.1, KI270803.1, KI270892.1, KI270897.1, KQ458383.1, KQ458385.1, ML143355.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, ML143352.1, ML143359.1, GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270764.1, KI270813.1, KI270816.1, KI270819.1, KI270832.1, KI270844.1, KI270850.1, KI270894.1, KN196484.1, KN196487.1, KN538364.1, ML143366.1, ML143382.1, GL000216.2, KI270587.1, KI270720.1, KI270892.1, KI270897.1
    ##   - in 'y': GL339449.2, GL383526.1, GL383581.2, KI270731.1, KI270743.1, KI270848.1, KI270856.1, KI270861.1, KN538370.1, KV880764.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, ML143352.1, ML143359.1, GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270764.1, KI270813.1, KI270816.1, KI270819.1, KI270832.1, KI270844.1, KI270850.1, KI270894.1, KN196484.1, KN196487.1, KN538364.1, ML143366.1, ML143382.1, GL000216.2, KI270587.1, KI270720.1, KI270892.1, KI270897.1
    ##   - in 'y': GL339449.2, GL383526.1, GL383581.2, KI270731.1, KI270743.1, KI270848.1, KI270856.1, KI270861.1, KN538370.1, KV880764.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270712.1, ML143352.1, ML143359.1, GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270707.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270754.1, KI270764.1, KI270816.1, KI270844.1, KI270894.1, KN196484.1, KN196487.1, ML143353.1, ML143382.1, GL000216.2, KI270587.1, KI270720.1, KI270892.1, KI270897.1, KQ458385.1, ML143355.1, GL383526.1, GL383581.2, KI270731.1, KI270743.1, KI270856.1, KI270861.1, KN538370.1, ML143379.1
    ##   - in 'y': GL000254.2, GL383574.1, GL877875.1, KI270310.1, KI270784.1, KI270792.1, KI270876.1, KV880768.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270712.1, ML143352.1, ML143359.1, GL000194.1, GL000218.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270707.1, KI270709.1, KI270726.1, KI270729.1, KI270733.1, KI270751.1, KI270754.1, KI270764.1, KI270816.1, KI270844.1, KI270894.1, KN196484.1, KN196487.1, ML143353.1, ML143382.1, GL000216.2, KI270587.1, KI270720.1, KI270892.1, KI270897.1, KQ458385.1, ML143355.1, GL383526.1, GL383581.2, KI270731.1, KI270743.1, KI270856.1, KI270861.1, KN538370.1, ML143379.1
    ##   - in 'y': GL000254.2, GL383574.1, GL877875.1, KI270310.1, KI270784.1, KI270792.1, KI270876.1, KV880768.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000256.2, GL339449.2, GL383578.2, GL949752.1, KI270706.1, KI270728.1, KI270733.1, KI270744.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KZ559103.1, ML143377.1, ML143380.1
    ##   - in 'y': GL000214.1, KI270836.1, KI270849.1, KI270879.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000253.2, GL000256.2, GL339449.2, GL383578.2, GL949752.1, KI270706.1, KI270728.1, KI270733.1, KI270744.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KZ559103.1, ML143377.1, ML143380.1
    ##   - in 'y': GL000214.1, KI270836.1, KI270849.1, KI270879.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, GL383578.2, KI270728.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KZ559103.1, KI270836.1, KI270849.1, KI270879.1
    ##   - in 'y': GL000009.2, GL000220.1, GL000225.1, KI270310.1, KI270438.1, KI270709.1, KQ458385.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, GL383578.2, KI270728.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KZ559103.1, KI270836.1, KI270849.1, KI270879.1
    ##   - in 'y': GL000009.2, GL000220.1, GL000225.1, KI270310.1, KI270438.1, KI270709.1, KQ458385.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL339449.2, GL383578.2, GL949752.1, KI270706.1, KI270733.1, KI270744.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KV766198.1, KZ559103.1, ML143377.1, ML143380.1, GL000214.1, KI270836.1, KI270849.1, KI270879.1, GL000009.2, GL000220.1, GL000225.1, KI270310.1, KI270438.1, KI270709.1, KQ458385.1
    ##   - in 'y': KZ559112.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000256.2, GL339449.2, GL383578.2, GL949752.1, KI270706.1, KI270733.1, KI270744.1, KI270754.1, KI270765.1, KI270784.1, KI270857.1, KI270899.1, KI270905.1, KV766198.1, KZ559103.1, ML143377.1, ML143380.1, GL000214.1, KI270836.1, KI270849.1, KI270879.1, GL000009.2, GL000220.1, GL000225.1, KI270310.1, KI270438.1, KI270709.1, KQ458385.1
    ##   - in 'y': KZ559112.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000252.2, GL339449.2, KI270466.1, KI270467.1, KI270712.1, KI270765.1, KI270849.1, KI270855.1, KI270871.1, KI270878.1, KN196479.1, KN538370.1, KQ031389.1, KZ208912.1, KZ559109.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000220.1, GL000225.1, GL000257.2, JH159146.1, KI270310.1, KI270438.1, KI270442.1, KI270707.1, KI270709.1, KI270719.1, KI270726.1, KI270730.1, KI270734.1, KI270748.1, KI270751.1, KI270757.1, KI270762.1, KI270792.1, KI270813.1, KI270844.1, KI270846.1, KI270848.1, KI270856.1, KI270894.1, KI270897.1, KI270904.1, KI270925.1, KN196482.1, KN538372.1, KQ458385.1, KZ559105.1, KZ559112.1, ML143344.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000252.2, GL339449.2, KI270466.1, KI270467.1, KI270712.1, KI270765.1, KI270849.1, KI270855.1, KI270871.1, KI270878.1, KN196479.1, KN538370.1, KQ031389.1, KZ208912.1, KZ559109.1, ML143352.1, ML143353.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000220.1, GL000225.1, GL000257.2, JH159146.1, KI270310.1, KI270438.1, KI270442.1, KI270707.1, KI270709.1, KI270719.1, KI270726.1, KI270730.1, KI270734.1, KI270748.1, KI270751.1, KI270757.1, KI270762.1, KI270792.1, KI270813.1, KI270844.1, KI270846.1, KI270848.1, KI270856.1, KI270894.1, KI270897.1, KI270904.1, KI270925.1, KN196482.1, KN538372.1, KQ458385.1, KZ559105.1, KZ559112.1, ML143344.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270538.1, ML143366.1
    ##   - in 'y': GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383579.2, GL877875.1, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270713.1, KI270714.1, KI270728.1, KI270745.1, KI270765.1, KI270787.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270905.1, KN538364.1, KQ031389.1, KV575244.1, KV766193.1, KZ208907.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143353.1, ML143355.1, ML143369.1, ML143371.1, ML143372.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270538.1, ML143366.1
    ##   - in 'y': GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383579.2, GL877875.1, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270713.1, KI270714.1, KI270728.1, KI270745.1, KI270765.1, KI270787.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270905.1, KN538364.1, KQ031389.1, KV575244.1, KV766193.1, KZ208907.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143353.1, ML143355.1, ML143369.1, ML143371.1, ML143372.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270744.1, KI270879.1, KI270908.1, KV766198.1, KV880768.1, ML143366.1, GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383579.2, GL877875.1, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270713.1, KI270714.1, KI270728.1, KI270745.1, KI270765.1, KI270787.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270905.1, KN538364.1, KQ031389.1, KV575244.1, KV766193.1, KZ208907.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143353.1, ML143369.1, ML143371.1, ML143372.1, ML143379.1
    ##   - in 'y': KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270744.1, KI270879.1, KI270908.1, KV766198.1, KV880768.1, ML143366.1, GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383556.1, GL383579.2, GL877875.1, GL949752.1, KI270442.1, KI270706.1, KI270707.1, KI270713.1, KI270714.1, KI270728.1, KI270745.1, KI270765.1, KI270787.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270905.1, KN538364.1, KQ031389.1, KV575244.1, KV766193.1, KZ208907.1, KZ208915.1, KZ559105.1, KZ559109.1, KZ559112.1, ML143353.1, ML143369.1, ML143371.1, ML143372.1, ML143379.1
    ##   - in 'y': KI270466.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270538.1, KI270908.1, KV880768.1, ML143366.1, GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000254.2, GL383556.1, GL383579.2, GL949752.1, KI270707.1, KI270714.1, KI270830.1, KI270849.1, KI270857.1, KQ031389.1, KV575244.1, KZ559105.1, ML143353.1, ML143372.1, ML143379.1, KI270466.1
    ##   - in 'y': GL383519.1, GL949747.2, KI270783.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270538.1, KI270908.1, KV880768.1, ML143366.1, GL000009.2, GL000219.1, GL000221.1, GL000251.2, GL000254.2, GL383556.1, GL383579.2, GL949752.1, KI270707.1, KI270714.1, KI270830.1, KI270849.1, KI270857.1, KQ031389.1, KV575244.1, KZ559105.1, ML143353.1, ML143372.1, ML143379.1, KI270466.1
    ##   - in 'y': GL383519.1, GL949747.2, KI270783.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL949753.2, KI270330.1, KI270337.1, KI270467.1, KI270712.1, KI270871.1, KI270872.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143371.1, ML143379.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000225.1, KI270706.1, KI270713.1, KI270744.1, KI270754.1, KI270784.1, KI270845.1, KI270853.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL949753.2, KI270330.1, KI270337.1, KI270467.1, KI270712.1, KI270871.1, KI270872.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143371.1, ML143379.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000225.1, KI270706.1, KI270713.1, KI270744.1, KI270754.1, KI270784.1, KI270845.1, KI270853.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270718.1, KI270728.1, KI270765.1, KI270782.1, KV880764.1, ML143353.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL949746.1, GL949752.1, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270726.1, KI270733.1, KI270736.1, KI270742.1, KI270745.1, KI270754.1, KI270770.1, KI270784.1, KI270809.1, KI270816.1, KI270819.1, KI270830.1, KI270832.1, KI270849.1, KI270851.1, KI270856.1, KI270860.1, KI270861.1, KI270868.1, KI270879.1, KI270897.1, KI270908.1, KN196479.1, KN196487.1, KN538364.1, KV766198.1, KV880768.1, KZ208914.1, KZ208915.1, KZ208921.1, KZ559112.1, ML143345.1, ML143370.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270718.1, KI270728.1, KI270765.1, KI270782.1, KV880764.1, ML143353.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL000251.2, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL339449.2, GL383563.3, GL949746.1, GL949752.1, KI270706.1, KI270711.1, KI270713.1, KI270714.1, KI270726.1, KI270733.1, KI270736.1, KI270742.1, KI270745.1, KI270754.1, KI270770.1, KI270784.1, KI270809.1, KI270816.1, KI270819.1, KI270830.1, KI270832.1, KI270849.1, KI270851.1, KI270856.1, KI270860.1, KI270861.1, KI270868.1, KI270879.1, KI270897.1, KI270908.1, KN196479.1, KN196487.1, KN538364.1, KV766198.1, KV880768.1, KZ208914.1, KZ208915.1, KZ208921.1, KZ559112.1, ML143345.1, ML143370.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270745.1, KI270908.1, KI270937.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL000224.1, GL000255.2, GL949752.1, KI270712.1, KI270721.1, KI270731.1, KI270738.1, KI270742.1, KI270744.1, KI270754.1, KI270770.1, KI270830.1, KI270850.1, KI270879.1, KN196487.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270745.1, KI270908.1, KI270937.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000220.1, GL000224.1, GL000255.2, GL949752.1, KI270712.1, KI270721.1, KI270731.1, KI270738.1, KI270742.1, KI270744.1, KI270754.1, KI270770.1, KI270830.1, KI270850.1, KI270879.1, KN196487.1, ML143371.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383526.1, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270908.1, ML143350.1, ML143352.1, ML143375.1, chr11, chr16, chr18, chr19, chr20, chr21, chr22, chr3, chr5, chr7, chr8, chr9, chrX
    ##   - in 'y': GL000224.1, GL000225.1, KI270438.1, KN196487.1, KZ559100.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383526.1, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270908.1, ML143350.1, ML143352.1, ML143375.1, chr11, chr16, chr18, chr19, chr20, chr21, chr22, chr3, chr5, chr7, chr8, chr9, chrX
    ##   - in 'y': GL000224.1, GL000225.1, KI270438.1, KN196487.1, KZ559100.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000251.2, KI270467.1, KI270719.1, KI270721.1, KI270878.1, ML143369.1
    ##   - in 'y': GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270830.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL000251.2, KI270467.1, KI270719.1, KI270721.1, KI270878.1, ML143369.1
    ##   - in 'y': GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270830.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL383527.1, KI270442.1, KI270467.1, KI270709.1, KI270719.1, KI270721.1, KI270783.1, KI270803.1, KI270861.1, KV880768.1, KZ559103.1, GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1
    ##   - in 'y': GL383578.2, KI270765.1, KI270784.1, KI270832.1, KI270849.1, KI270850.1, KI270857.1, KZ208912.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL383527.1, KI270442.1, KI270467.1, KI270709.1, KI270719.1, KI270721.1, KI270783.1, KI270803.1, KI270861.1, KV880768.1, KZ559103.1, GL000009.2, GL000205.2, GL000216.2, GL000220.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1
    ##   - in 'y': GL383578.2, KI270765.1, KI270784.1, KI270832.1, KI270849.1, KI270850.1, KI270857.1, KZ208912.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL383527.1, GL877875.1, KI270467.1, KI270721.1, KI270728.1, KI270783.1, KI270803.1, KI270831.1, KI270861.1, KI270878.1, KI270908.1, KV766198.1, KZ559103.1, KZ559112.1, ML143366.1, GL000009.2, GL000205.2, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1, GL383578.2, KI270765.1, KI270784.1, KI270849.1, KI270850.1, KZ208912.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': KI270813.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000218.1, GL000221.1, GL000225.1, GL000251.2, GL383527.1, GL877875.1, KI270467.1, KI270721.1, KI270728.1, KI270783.1, KI270803.1, KI270831.1, KI270861.1, KI270878.1, KI270908.1, KV766198.1, KZ559103.1, KZ559112.1, ML143366.1, GL000009.2, GL000205.2, GL000224.1, GL383522.1, GL383567.1, KI270311.1, KI270330.1, KI270507.1, KI270516.1, KI270538.1, KI270590.1, KI270707.1, KI270712.1, KI270714.1, KI270723.1, KI270729.1, KI270732.1, KI270733.1, KI270738.1, KI270746.1, KI270751.1, KI270754.1, KI270757.1, KI270767.1, KI270787.1, KI270816.1, KI270836.1, KI270851.1, KI270905.1, KN196472.1, KN196487.1, KN538372.1, KQ983257.1, KV766196.1, ML143371.1, GL383578.2, KI270765.1, KI270784.1, KI270849.1, KI270850.1, KZ208912.1, ML143352.1, ML143355.1, ML143375.1, ML143379.1
    ##   - in 'y': KI270813.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270726.1, KI270731.1, KI270799.1, KI270831.1, KI270844.1, KI270856.1, KI270868.1, KI270876.1, KI270903.1, KN196472.1, KQ458384.1, KV880764.1, KZ208912.1, ML143345.1, ML143355.1, ML143377.1
    ##   - in 'y': GL000224.1, GL383556.1, KI270310.1, KI270816.1, KI270849.1, KN538361.1, KQ031389.1, KZ208913.1, KZ208922.1, KZ559103.1, KZ559109.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270726.1, KI270731.1, KI270799.1, KI270831.1, KI270844.1, KI270856.1, KI270868.1, KI270876.1, KI270903.1, KN196472.1, KQ458384.1, KV880764.1, KZ208912.1, ML143345.1, ML143355.1, ML143377.1
    ##   - in 'y': GL000224.1, GL383556.1, KI270310.1, KI270816.1, KI270849.1, KN538361.1, KQ031389.1, KZ208913.1, KZ208922.1, KZ559103.1, KZ559109.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, ML143380.1
    ##   - in 'y': GL000253.2, KI270714.1, KV766198.1, chr10, chr21
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, ML143380.1
    ##   - in 'y': GL000253.2, KI270714.1, KV766198.1, chr10, chr21
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, ML143380.1, GL000253.2, KI270714.1, KV766198.1
    ##   - in 'y': GL000224.1, KI270892.1, KI270907.1, KN538364.1, ML143375.1, chr18
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, ML143380.1, GL000253.2, KI270714.1, KV766198.1
    ##   - in 'y': GL000224.1, KI270892.1, KI270907.1, KN538364.1, ML143375.1, chr18
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, ML143380.1, KI270714.1, KV766198.1, chr21, GL000224.1, KI270892.1, KI270907.1, KN538364.1, ML143375.1
    ##   - in 'y': GL000256.2, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, ML143380.1, KI270714.1, KV766198.1, chr21, GL000224.1, KI270892.1, KI270907.1, KN538364.1, ML143375.1
    ##   - in 'y': GL000256.2, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270844.1
    ##   - in 'y': GL000216.2, GL000252.2, GL383563.3, KI270731.1, KI270851.1, KN196487.1, KZ208914.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270844.1
    ##   - in 'y': GL000216.2, GL000252.2, GL383563.3, KI270731.1, KI270851.1, KN196487.1, KZ208914.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383545.1, GL383579.2, GL949746.1, KI270734.1, KI270745.1, KI270762.1, KI270780.1, KI270795.1, KI270809.1, KI270830.1, KI270831.1, KI270832.1, KI270844.1, KI270861.1, KI270904.1, KQ090026.1, KZ208913.1, KZ559105.1, GL000252.2, GL383563.3, KI270731.1, KI270851.1
    ##   - in 'y': GL000224.1, GL877875.1, KI270320.1, KI270709.1, KI270746.1, KI270751.1, KI270770.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383545.1, GL383579.2, GL949746.1, KI270734.1, KI270745.1, KI270762.1, KI270780.1, KI270795.1, KI270809.1, KI270830.1, KI270831.1, KI270832.1, KI270844.1, KI270861.1, KI270904.1, KQ090026.1, KZ208913.1, KZ559105.1, GL000252.2, GL383563.3, KI270731.1, KI270851.1
    ##   - in 'y': GL000224.1, GL877875.1, KI270320.1, KI270709.1, KI270746.1, KI270751.1, KI270770.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL383579.2, KI270438.1, KI270723.1, KI270762.1, KI270780.1, KI270809.1, KI270844.1, KI270904.1, KZ559105.1, GL000216.2, GL000252.2, GL383563.3, KI270731.1, KI270851.1, KN196487.1, ML143366.1, GL877875.1, KI270320.1, KI270746.1, KI270751.1
    ##   - in 'y': KI270849.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL383579.2, KI270438.1, KI270723.1, KI270762.1, KI270780.1, KI270809.1, KI270844.1, KI270904.1, KZ559105.1, GL000216.2, GL000252.2, GL383563.3, KI270731.1, KI270851.1, KN196487.1, ML143366.1, GL877875.1, KI270320.1, KI270746.1, KI270751.1
    ##   - in 'y': KI270849.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270908.1, ML143352.1
    ##   - in 'y': GL000218.1, GL000220.1, GL000224.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270908.1, ML143352.1
    ##   - in 'y': GL000218.1, GL000220.1, GL000224.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270733.1, KI270754.1, KI270765.1, KI270782.1, KI270813.1, ML143352.1, ML143355.1, ML143372.1
    ##   - in 'y': GL000253.2, GL000255.2, KI270706.1, KI270742.1, KI270744.1, KI270880.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270733.1, KI270754.1, KI270765.1, KI270782.1, KI270813.1, ML143352.1, ML143355.1, ML143372.1
    ##   - in 'y': GL000253.2, GL000255.2, KI270706.1, KI270742.1, KI270744.1, KI270880.1, KQ090026.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270330.1, KI270435.1, KI270712.1, KI270718.1, KI270729.1, KI270735.1, KI270736.1, KI270738.1, KI270765.1, KV766198.1, ML143352.1
    ##   - in 'y': GL000205.2, GL000221.1, GL000255.2, GL000256.2, GL949747.2, GL949752.1, KI270587.1, KI270714.1, KI270726.1, KI270743.1, KI270751.1, KI270754.1, KI270762.1, KI270767.1, KI270805.1, KI270816.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270846.1, KI270850.1, KI270851.1, KI270853.1, KI270861.1, KI270892.1, KI270897.1, KI270905.1, KI270908.1, KI270924.1, KN538361.1, KZ208915.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270330.1, KI270435.1, KI270712.1, KI270718.1, KI270729.1, KI270735.1, KI270736.1, KI270738.1, KI270765.1, KV766198.1, ML143352.1
    ##   - in 'y': GL000205.2, GL000221.1, GL000255.2, GL000256.2, GL949747.2, GL949752.1, KI270587.1, KI270714.1, KI270726.1, KI270743.1, KI270751.1, KI270754.1, KI270762.1, KI270767.1, KI270805.1, KI270816.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270846.1, KI270850.1, KI270851.1, KI270853.1, KI270861.1, KI270892.1, KI270897.1, KI270905.1, KI270908.1, KI270924.1, KN538361.1, KZ208915.1, ML143371.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000253.2, GL000256.2, KI270311.1, KI270330.1, KI270731.1, KI270738.1, KI270742.1, KI270744.1, KI270754.1, KI270816.1, KI270844.1, KI270856.1, KI270879.1, KI270880.1, KI270905.1, KI270908.1, KQ090026.1, KV575244.1, KV880763.1, KV880768.1, ML143345.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000253.2, GL000256.2, KI270311.1, KI270330.1, KI270731.1, KI270738.1, KI270742.1, KI270744.1, KI270754.1, KI270816.1, KI270844.1, KI270856.1, KI270879.1, KI270880.1, KI270905.1, KI270908.1, KQ090026.1, KV575244.1, KV880763.1, KV880768.1, ML143345.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383522.1, GL383545.1, GL383557.1, KI270707.1, KI270726.1, KI270734.1, KI270754.1, KI270779.1, KI270783.1, KI270787.1, KI270827.1, KI270831.1, KI270848.1, KI270861.1, KI270875.1, KI270892.1, KI270907.1, KN538361.1, KQ031389.1, KV880764.1, KZ208907.1, KZ559105.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000221.1, GL383582.2, KI270719.1, KI270804.1, KI270821.1, KI270822.1, KI270836.1, KI270849.1, KN196484.1, KQ458383.1, KQ458385.1, KQ983255.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383522.1, GL383545.1, GL383557.1, KI270707.1, KI270726.1, KI270734.1, KI270754.1, KI270779.1, KI270783.1, KI270787.1, KI270827.1, KI270831.1, KI270848.1, KI270861.1, KI270875.1, KI270892.1, KI270907.1, KN538361.1, KQ031389.1, KV880764.1, KZ208907.1, KZ559105.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000221.1, GL383582.2, KI270719.1, KI270804.1, KI270821.1, KI270822.1, KI270836.1, KI270849.1, KN196484.1, KQ458383.1, KQ458385.1, KQ983255.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270734.1, KN538370.1, KQ090027.1, KV880763.1, ML143359.1, ML143367.1
    ##   - in 'y': GL383557.1, GL877875.1, KI270804.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270734.1, KN538370.1, KQ090027.1, KV880763.1, ML143359.1, ML143367.1
    ##   - in 'y': GL383557.1, GL877875.1, KI270804.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383519.1, GL383522.1, GL383579.2, GL949746.1, KI270442.1, KI270720.1, KI270734.1, KI270750.1, KI270831.1, KI270844.1, KI270855.1, KI270860.1, KI270861.1, KI270866.1, KI270899.1, KN196479.1, KN538370.1, KQ090027.1, KQ458383.1, KV880763.1, KZ208914.1, KZ559105.1, ML143359.1, ML143367.1, ML143371.1, GL383557.1, GL877875.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1
    ##   - in 'y': KI270871.1, KI270936.1, KQ090026.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL383519.1, GL383522.1, GL383579.2, GL949746.1, KI270442.1, KI270720.1, KI270734.1, KI270750.1, KI270831.1, KI270844.1, KI270855.1, KI270860.1, KI270861.1, KI270866.1, KI270899.1, KN196479.1, KN538370.1, KQ090027.1, KQ458383.1, KV880763.1, KZ208914.1, KZ559105.1, ML143359.1, ML143367.1, ML143371.1, GL383557.1, GL877875.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1
    ##   - in 'y': KI270871.1, KI270936.1, KQ090026.1, ML143345.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000221.1, GL000254.2, GL383519.1, GL383522.1, GL949746.1, KI270442.1, KI270706.1, KI270720.1, KI270721.1, KI270734.1, KI270750.1, KI270782.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270844.1, KI270850.1, KI270855.1, KI270860.1, KI270861.1, KI270866.1, KI270868.1, KI270897.1, KI270899.1, KI270908.1, KN196479.1, KN538364.1, KN538370.1, KN538372.1, KQ090027.1, KQ458383.1, KV880763.1, KZ208914.1, KZ559105.1, ML143359.1, ML143367.1, GL383557.1, GL877875.1, KI270804.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1, KI270936.1, KQ090026.1
    ##   - in 'y': KI270719.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000221.1, GL000254.2, GL383519.1, GL383522.1, GL949746.1, KI270442.1, KI270706.1, KI270720.1, KI270721.1, KI270734.1, KI270750.1, KI270782.1, KI270816.1, KI270830.1, KI270831.1, KI270832.1, KI270844.1, KI270850.1, KI270855.1, KI270860.1, KI270861.1, KI270866.1, KI270868.1, KI270897.1, KI270899.1, KI270908.1, KN196479.1, KN538364.1, KN538370.1, KN538372.1, KQ090027.1, KQ458383.1, KV880763.1, KZ208914.1, KZ559105.1, ML143359.1, ML143367.1, GL383557.1, GL877875.1, KI270804.1, KI270809.1, KI270851.1, KI270865.1, KI270895.1, KN538373.1, KQ458384.1, KV575243.1, KV766194.1, KZ208907.1, KZ208917.1, ML143343.1, KI270936.1, KQ090026.1
    ##   - in 'y': KI270719.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270721.1, KI270733.1, KI270750.1, KI270765.1, KI270851.1, KI270855.1, KI270868.1, KI270895.1, KN538370.1, ML143344.1, ML143352.1
    ##   - in 'y': GL383557.1, KI270583.1, KI270707.1, KI270734.1, KI270743.1, KI270809.1, KI270821.1, KI270830.1, KI270844.1, KI270856.1, KI270897.1, KQ090026.1, KZ559103.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270712.1, KI270721.1, KI270733.1, KI270750.1, KI270765.1, KI270851.1, KI270855.1, KI270868.1, KI270895.1, KN538370.1, ML143344.1, ML143352.1
    ##   - in 'y': GL383557.1, KI270583.1, KI270707.1, KI270734.1, KI270743.1, KI270809.1, KI270821.1, KI270830.1, KI270844.1, KI270856.1, KI270897.1, KQ090026.1, KZ559103.1, ML143367.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270765.1, KI270803.1, KZ208913.1, KZ559109.1, ML143352.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000221.1, GL339449.2, GL383574.1, KI270330.1, KI270590.1, KI270706.1, KI270711.1, KI270719.1, KI270762.1, KI270813.1, KI270827.1, KQ031389.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270765.1, KI270803.1, KZ208913.1, KZ559109.1, ML143352.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000221.1, GL339449.2, GL383574.1, KI270330.1, KI270590.1, KI270706.1, KI270711.1, KI270719.1, KI270762.1, KI270813.1, KI270827.1, KQ031389.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1
    ##   - in 'y': GL000218.1, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1
    ##   - in 'y': GL000218.1, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000220.1, GL000225.1, KI270442.1, KI270706.1, KI270744.1, KZ559112.1, ML143372.1, GL000218.1, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   - in 'y': KI270842.1, KI270850.1, KI270908.1, KQ090026.1, ML143345.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000220.1, GL000225.1, KI270442.1, KI270706.1, KI270744.1, KZ559112.1, ML143372.1, GL000218.1, GL000256.2, KI270438.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1
    ##   - in 'y': KI270842.1, KI270850.1, KI270908.1, KQ090026.1, ML143345.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000220.1, GL000225.1, KI270442.1, KI270706.1, KZ559112.1, ML143372.1, GL000218.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1, KI270850.1, KI270908.1, KQ090026.1, ML143345.1
    ##   - in 'y': GL000221.1, GL000226.1, KI270467.1, KI270709.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL000220.1, GL000225.1, KI270442.1, KI270706.1, KZ559112.1, ML143372.1, GL000218.1, KI270714.1, KI270905.1, KQ983257.1, KZ208915.1, KI270850.1, KI270908.1, KQ090026.1, ML143345.1
    ##   - in 'y': GL000221.1, GL000226.1, KI270467.1, KI270709.1, KI270728.1, KI270729.1, KI270733.1, KI270736.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270709.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, GL000251.2, GL000255.2, GL383580.2, GL949752.1, KI270711.1, KI270721.1, KI270726.1, KI270742.1, KI270757.1, KI270780.1, KI270792.1, KI270816.1, KI270821.1, KI270844.1, KI270850.1, KI270853.1, KI270856.1, KI270879.1, KI270880.1, KI270903.1, KI270905.1, KN196484.1, KN196487.1, KQ458383.1, KQ983257.1, KV575244.1, KZ208917.1, KZ559103.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270709.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000224.1, GL000225.1, GL000251.2, GL000255.2, GL383580.2, GL949752.1, KI270711.1, KI270721.1, KI270726.1, KI270742.1, KI270757.1, KI270780.1, KI270792.1, KI270816.1, KI270821.1, KI270844.1, KI270850.1, KI270853.1, KI270856.1, KI270879.1, KI270880.1, KI270903.1, KI270905.1, KN196484.1, KN196487.1, KQ458383.1, KQ983257.1, KV575244.1, KZ208917.1, KZ559103.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000216.2, KI270438.1, KI270538.1, KI270709.1
    ##   - in 'y': GL000220.1, GL000250.2, GL000251.2, GL000253.2, GL000257.2, GL339449.2, GL383527.1, GL383545.1, GL383556.1, GL383563.3, GL383580.2, GL877875.1, GL949752.1, KI270442.1, KI270708.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270722.1, KI270725.1, KI270728.1, KI270733.1, KI270734.1, KI270761.1, KI270770.1, KI270773.1, KI270780.1, KI270783.1, KI270792.1, KI270802.1, KI270803.1, KI270804.1, KI270813.1, KI270819.1, KI270821.1, KI270830.1, KI270849.1, KI270850.1, KI270851.1, KI270861.1, KI270862.1, KI270880.1, KI270892.1, KI270894.1, KI270896.1, KI270907.1, KI270938.1, KN196484.1, KQ458383.1, KQ458384.1, KV766198.1, KV880768.1, KZ208913.1, KZ208921.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000214.1, GL000216.2, KI270438.1, KI270538.1, KI270709.1
    ##   - in 'y': GL000220.1, GL000250.2, GL000251.2, GL000253.2, GL000257.2, GL339449.2, GL383527.1, GL383545.1, GL383556.1, GL383563.3, GL383580.2, GL877875.1, GL949752.1, KI270442.1, KI270708.1, KI270712.1, KI270713.1, KI270714.1, KI270717.1, KI270719.1, KI270722.1, KI270725.1, KI270728.1, KI270733.1, KI270734.1, KI270761.1, KI270770.1, KI270773.1, KI270780.1, KI270783.1, KI270792.1, KI270802.1, KI270803.1, KI270804.1, KI270813.1, KI270819.1, KI270821.1, KI270830.1, KI270849.1, KI270850.1, KI270851.1, KI270861.1, KI270862.1, KI270880.1, KI270892.1, KI270894.1, KI270896.1, KI270907.1, KI270938.1, KN196484.1, KQ458383.1, KQ458384.1, KV766198.1, KV880768.1, KZ208913.1, KZ208921.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000257.2, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270746.1, KI270857.1, KI270908.1, KN196472.1, KV766198.1, KZ208912.1, ML143380.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000224.1, GL000225.1, KI270330.1, KI270538.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270743.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000257.2, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270746.1, KI270857.1, KI270908.1, KN196472.1, KV766198.1, KZ208912.1, ML143380.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000224.1, GL000225.1, KI270330.1, KI270538.1, KI270731.1, KI270733.1, KI270734.1, KI270736.1, KI270743.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000257.2, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270746.1, KI270908.1, KN196472.1, KV766198.1, KZ208912.1, GL000214.1, GL000224.1, KI270330.1, KI270538.1, KI270731.1, KI270736.1, KI270743.1
    ##   - in 'y': GL339449.2, GL383555.2, GL383580.2, KI270787.1, KI270861.1, KI270879.1, KI270905.1, KN538364.1, KQ090027.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000257.2, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270722.1, KI270746.1, KI270908.1, KN196472.1, KV766198.1, KZ208912.1, GL000214.1, GL000224.1, KI270330.1, KI270538.1, KI270731.1, KI270736.1, KI270743.1
    ##   - in 'y': GL339449.2, GL383555.2, GL383580.2, KI270787.1, KI270861.1, KI270879.1, KI270905.1, KN538364.1, KQ090027.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL582966.2, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270722.1, KI270746.1, KI270908.1, KN196487.1, GL000214.1, GL000224.1, GL000225.1, KI270330.1, KI270736.1, KI270787.1, KN538364.1
    ##   - in 'y': GL000226.1, GL383563.3, KI270714.1, KI270808.1, KI270871.1, KQ031384.1, KV575244.1, KV880764.1, KZ559112.1, ML143345.1, ML143375.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL582966.2, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270722.1, KI270746.1, KI270908.1, KN196487.1, GL000214.1, GL000224.1, GL000225.1, KI270330.1, KI270736.1, KI270787.1, KN538364.1
    ##   - in 'y': GL000226.1, GL383563.3, KI270714.1, KI270808.1, KI270871.1, KQ031384.1, KV575244.1, KV880764.1, KZ559112.1, ML143345.1, ML143375.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, KI270765.1, KI270830.1, KI270842.1, KI270905.1, KN538370.1, KV766198.1, KV880764.1, ML143352.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000225.1, KI270310.1, KI270465.1, KI270728.1, KI270729.1, KI270733.1, KI270744.1, KI270754.1, KI270850.1, KN196487.1, KN538364.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL339449.2, KI270765.1, KI270830.1, KI270842.1, KI270905.1, KN538370.1, KV766198.1, KV880764.1, ML143352.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000224.1, GL000225.1, KI270310.1, KI270465.1, KI270728.1, KI270729.1, KI270733.1, KI270744.1, KI270754.1, KI270850.1, KN196487.1, KN538364.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, JH636055.2, KI270728.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000205.2, GL000214.1, GL000220.1, GL000224.1, GL000252.2, GL000256.2, GL383556.1, KI270330.1, KI270442.1, KI270742.1, KI270751.1, KI270754.1, KI270759.1, KI270763.1, KI270782.1, KI270802.1, KI270808.1, KI270835.1, KI270851.1, KI270878.1, KI270897.1, KI270908.1, KI270926.1, KI270936.1, KI270937.1, KQ031389.1, KQ090017.1, KQ458385.1, KZ208913.1, KZ559103.1, KZ559105.1, ML143366.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, JH636055.2, KI270728.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000205.2, GL000214.1, GL000220.1, GL000224.1, GL000252.2, GL000256.2, GL383556.1, KI270330.1, KI270442.1, KI270742.1, KI270751.1, KI270754.1, KI270759.1, KI270763.1, KI270782.1, KI270802.1, KI270808.1, KI270835.1, KI270851.1, KI270878.1, KI270897.1, KI270908.1, KI270926.1, KI270936.1, KI270937.1, KQ031389.1, KQ090017.1, KQ458385.1, KZ208913.1, KZ559103.1, KZ559105.1, ML143366.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, JH636055.2, ML143352.1, GL000008.2, GL383556.1, KI270759.1, KI270835.1, KQ090017.1, KZ559103.1
    ##   - in 'y': GL000194.1, GL000216.2, GL000218.1, GL000225.1, GL000253.2, GL000254.2, GL383563.3, GL383567.1, GL383580.2, JH159146.1, KI270713.1, KI270714.1, KI270723.1, KI270734.1, KI270736.1, KI270745.1, KI270792.1, KI270804.1, KI270809.1, KI270816.1, KI270818.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270836.1, KI270850.1, KI270856.1, KI270861.1, KI270868.1, KI270873.1, KI270875.1, KI270880.1, KI270892.1, KI270901.1, KI270903.1, KN538361.1, KN538364.1, KN538372.1, KQ090026.1, KV575244.1, KV766192.1, KV880768.1, KZ208922.1, ML143345.1, ML143353.1, ML143367.1, ML143375.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000221.1, GL383522.1, JH636055.2, ML143352.1, GL000008.2, GL383556.1, KI270759.1, KI270835.1, KQ090017.1, KZ559103.1
    ##   - in 'y': GL000194.1, GL000216.2, GL000218.1, GL000225.1, GL000253.2, GL000254.2, GL383563.3, GL383567.1, GL383580.2, JH159146.1, KI270713.1, KI270714.1, KI270723.1, KI270734.1, KI270736.1, KI270745.1, KI270792.1, KI270804.1, KI270809.1, KI270816.1, KI270818.1, KI270819.1, KI270821.1, KI270827.1, KI270830.1, KI270831.1, KI270832.1, KI270836.1, KI270850.1, KI270856.1, KI270861.1, KI270868.1, KI270873.1, KI270875.1, KI270880.1, KI270892.1, KI270901.1, KI270903.1, KN538361.1, KN538364.1, KN538372.1, KQ090026.1, KV575244.1, KV766192.1, KV880768.1, KZ208922.1, ML143345.1, ML143353.1, ML143367.1, ML143375.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383545.1, JH636055.2, KI270709.1, KI270712.1, KI270717.1, KI270719.1, KI270733.1, KI270765.1, KI270772.1, KI270842.1, KI270894.1, GL000205.2, GL000220.1, GL000252.2, GL383556.1, KI270330.1, KI270751.1, KI270759.1, KI270763.1, KI270802.1, KI270808.1, KI270835.1, KI270926.1, KI270936.1, KI270937.1, KQ090017.1, KZ208913.1, GL000216.2, GL000225.1, GL383563.3, GL383567.1, GL383580.2, JH159146.1, KI270723.1, KI270736.1, KI270804.1, KI270818.1, KI270827.1, KI270836.1, KI270873.1, KI270875.1, KI270892.1, KI270901.1, KN538372.1, KV766192.1, ML143367.1, ML143378.1
    ##   - in 'y': GL000251.2, GL000257.2, GL383520.2, GL383578.2, GL877875.1, GL949746.1, KI270706.1, KI270707.1, KI270711.1, KI270718.1, KI270720.1, KI270721.1, KI270726.1, KI270743.1, KI270762.1, KI270770.1, KI270779.1, KI270783.1, KI270784.1, KI270787.1, KI270844.1, KI270848.1, KI270849.1, KI270853.1, KI270869.1, KI270876.1, KI270879.1, KI270904.1, KN196472.1, KN538370.1, KQ458383.1, KV766193.1, KV766198.1, KZ208907.1, KZ208914.1, KZ559109.1, KZ559112.1, ML143360.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383545.1, JH636055.2, KI270709.1, KI270712.1, KI270717.1, KI270719.1, KI270733.1, KI270765.1, KI270772.1, KI270842.1, KI270894.1, GL000205.2, GL000220.1, GL000252.2, GL383556.1, KI270330.1, KI270751.1, KI270759.1, KI270763.1, KI270802.1, KI270808.1, KI270835.1, KI270926.1, KI270936.1, KI270937.1, KQ090017.1, KZ208913.1, GL000216.2, GL000225.1, GL383563.3, GL383567.1, GL383580.2, JH159146.1, KI270723.1, KI270736.1, KI270804.1, KI270818.1, KI270827.1, KI270836.1, KI270873.1, KI270875.1, KI270892.1, KI270901.1, KN538372.1, KV766192.1, ML143367.1, ML143378.1
    ##   - in 'y': GL000251.2, GL000257.2, GL383520.2, GL383578.2, GL877875.1, GL949746.1, KI270706.1, KI270707.1, KI270711.1, KI270718.1, KI270720.1, KI270721.1, KI270726.1, KI270743.1, KI270762.1, KI270770.1, KI270779.1, KI270783.1, KI270784.1, KI270787.1, KI270844.1, KI270848.1, KI270849.1, KI270853.1, KI270869.1, KI270876.1, KI270879.1, KI270904.1, KN196472.1, KN538370.1, KQ458383.1, KV766193.1, KV766198.1, KZ208907.1, KZ208914.1, KZ559109.1, KZ559112.1, ML143360.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270737.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000255.2, KI270590.1, KI270712.1, KI270729.1, KI270742.1, KI270746.1, KI270754.1, KI270757.1, KI270868.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270737.1
    ##   - in 'y': GL000194.1, GL000218.1, GL000255.2, KI270590.1, KI270712.1, KI270729.1, KI270742.1, KI270746.1, KI270754.1, KI270757.1, KI270868.1, KN538372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270709.1, KI270821.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, GL000225.1, GL000255.2, GL000256.2, GL949752.1, KI270712.1, KI270714.1, KI270728.1, KI270742.1, KI270744.1, KI270754.1, KI270830.1, KI270832.1, KI270853.1, KI270903.1, KI270905.1, KI270908.1, KN196484.1, KZ208917.1, KZ559103.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270709.1, KI270821.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000218.1, GL000219.1, GL000224.1, GL000225.1, GL000255.2, GL000256.2, GL949752.1, KI270712.1, KI270714.1, KI270728.1, KI270742.1, KI270744.1, KI270754.1, KI270830.1, KI270832.1, KI270853.1, KI270903.1, KI270905.1, KI270908.1, KN196484.1, KZ208917.1, KZ559103.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, GL000224.1, KI270538.1, KI270765.1, KI270907.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000221.1, GL000255.2, GL000256.2, GL383527.1, GL383578.2, GL949752.1, KI270320.1, KI270333.1, KI270337.1, KI270411.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270742.1, KI270745.1, KI270754.1, KI270762.1, KI270770.1, KI270784.1, KI270816.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270868.1, KI270879.1, KI270880.1, KI270892.1, KI270903.1, KN196484.1, KQ090026.1, KV766193.1, KV766198.1, KZ559105.1, ML143371.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000216.2, GL000224.1, KI270538.1, KI270765.1, KI270907.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000221.1, GL000255.2, GL000256.2, GL383527.1, GL383578.2, GL949752.1, KI270320.1, KI270333.1, KI270337.1, KI270411.1, KI270442.1, KI270466.1, KI270467.1, KI270706.1, KI270711.1, KI270714.1, KI270742.1, KI270745.1, KI270754.1, KI270762.1, KI270770.1, KI270784.1, KI270816.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270857.1, KI270868.1, KI270879.1, KI270880.1, KI270892.1, KI270903.1, KN196484.1, KQ090026.1, KV766193.1, KV766198.1, KZ559105.1, ML143371.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, GL000225.1, GL000250.2, GL000254.2, GL383545.1, GL877875.1, KI270712.1, KI270717.1, KI270718.1, KI270719.1, KI270720.1, KI270731.1, KI270751.1, KI270763.1, KI270772.1, KI270782.1, KI270792.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270827.1, KI270835.1, KI270842.1, KI270860.1, KI270873.1, KI270875.1, KI270876.1, KI270894.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538370.1, KN538372.1, KV880763.1, KV880764.1, KV880768.1, KZ208913.1, KZ208922.1, KZ559105.1, ML143352.1, ML143366.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000221.1, GL383520.2, GL383578.2, KI270438.1, KI270707.1, KI270711.1, KI270714.1, KI270721.1, KI270729.1, KI270736.1, KI270743.1, KI270770.1, KI270779.1, KI270783.1, KI270803.1, KI270819.1, KI270844.1, KI270849.1, KI270861.1, KI270869.1, KI270878.1, KN196479.1, KN538361.1, KQ090026.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, GL000225.1, GL000250.2, GL000254.2, GL383545.1, GL877875.1, KI270712.1, KI270717.1, KI270718.1, KI270719.1, KI270720.1, KI270731.1, KI270751.1, KI270763.1, KI270772.1, KI270782.1, KI270792.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270827.1, KI270835.1, KI270842.1, KI270860.1, KI270873.1, KI270875.1, KI270876.1, KI270894.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538370.1, KN538372.1, KV880763.1, KV880764.1, KV880768.1, KZ208913.1, KZ208922.1, KZ559105.1, ML143352.1, ML143366.1, ML143375.1
    ##   - in 'y': GL000009.2, GL000221.1, GL383520.2, GL383578.2, KI270438.1, KI270707.1, KI270711.1, KI270714.1, KI270721.1, KI270729.1, KI270736.1, KI270743.1, KI270770.1, KI270779.1, KI270783.1, KI270803.1, KI270819.1, KI270844.1, KI270849.1, KI270861.1, KI270869.1, KI270878.1, KN196479.1, KN538361.1, KQ090026.1, KZ559112.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, GL000225.1, GL000250.2, GL000254.2, GL383545.1, GL877875.1, KI270712.1, KI270717.1, KI270719.1, KI270720.1, KI270731.1, KI270733.1, KI270751.1, KI270754.1, KI270763.1, KI270772.1, KI270782.1, KI270784.1, KI270792.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270821.1, KI270827.1, KI270831.1, KI270835.1, KI270842.1, KI270860.1, KI270868.1, KI270873.1, KI270875.1, KI270876.1, KI270880.1, KI270894.1, KI270897.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538372.1, KV880763.1, KV880768.1, KZ208922.1, ML143352.1, ML143365.1, ML143372.1, GL000221.1, GL383578.2, KI270438.1, KI270729.1, KI270736.1, KI270743.1, KI270770.1, KI270779.1, KI270869.1, KN196479.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383581.2, KI270734.1, KI270816.1, KI270899.1, KI270904.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000224.1, GL000225.1, GL000250.2, GL000254.2, GL383545.1, GL877875.1, KI270712.1, KI270717.1, KI270719.1, KI270720.1, KI270731.1, KI270733.1, KI270751.1, KI270754.1, KI270763.1, KI270772.1, KI270782.1, KI270784.1, KI270792.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270821.1, KI270827.1, KI270831.1, KI270835.1, KI270842.1, KI270860.1, KI270868.1, KI270873.1, KI270875.1, KI270876.1, KI270880.1, KI270894.1, KI270897.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538372.1, KV880763.1, KV880768.1, KZ208922.1, ML143352.1, ML143365.1, ML143372.1, GL000221.1, GL383578.2, KI270438.1, KI270729.1, KI270736.1, KI270743.1, KI270770.1, KI270779.1, KI270869.1, KN196479.1
    ##   - in 'y': GL000008.2, GL383522.1, GL383581.2, KI270734.1, KI270816.1, KI270899.1, KI270904.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000214.1, GL000224.1, GL000225.1, GL000250.2, GL383545.1, KI270712.1, KI270717.1, KI270718.1, KI270719.1, KI270720.1, KI270731.1, KI270733.1, KI270751.1, KI270754.1, KI270763.1, KI270772.1, KI270782.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270827.1, KI270835.1, KI270842.1, KI270860.1, KI270873.1, KI270875.1, KI270894.1, KI270897.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538370.1, KN538372.1, KV880763.1, KV880764.1, KZ208922.1, ML143352.1, ML143353.1, ML143365.1, ML143375.1, GL383578.2, KI270438.1, KI270729.1, KI270736.1, KI270743.1, KI270779.1, KI270819.1, KI270869.1, KN196479.1, GL383581.2, KI270816.1, KI270899.1, KI270904.1
    ##   - in 'y': GL383556.1, GL383563.3, JH159146.1, KI270787.1, KZ559103.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL000214.1, GL000224.1, GL000225.1, GL000250.2, GL383545.1, KI270712.1, KI270717.1, KI270718.1, KI270719.1, KI270720.1, KI270731.1, KI270733.1, KI270751.1, KI270754.1, KI270763.1, KI270772.1, KI270782.1, KI270802.1, KI270804.1, KI270808.1, KI270809.1, KI270810.1, KI270827.1, KI270835.1, KI270842.1, KI270860.1, KI270873.1, KI270875.1, KI270894.1, KI270897.1, KI270924.1, KI270934.1, KI270936.1, KI270937.1, KN538370.1, KN538372.1, KV880763.1, KV880764.1, KZ208922.1, ML143352.1, ML143353.1, ML143365.1, ML143375.1, GL383578.2, KI270438.1, KI270729.1, KI270736.1, KI270743.1, KI270779.1, KI270819.1, KI270869.1, KN196479.1, GL383581.2, KI270816.1, KI270899.1, KI270904.1
    ##   - in 'y': GL383556.1, GL383563.3, JH159146.1, KI270787.1, KZ559103.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000224.1, GL000251.2, GL000255.2, GL383522.1, GL383556.1, GL383557.1, GL383578.2, GL877875.1, GL949752.1, KI270467.1, KI270706.1, KI270712.1, KI270714.1, KI270735.1, KI270745.1, KI270751.1, KI270765.1, KI270782.1, KI270804.1, KI270809.1, KI270845.1, KI270849.1, KI270850.1, KI270856.1, KI270861.1, KI270868.1, KI270878.1, KI270879.1, KI270897.1, KI270908.1, KI270934.1, KN196479.1, KV575244.1, KV766198.1, KV880764.1, KZ208912.1, KZ208915.1, ML143344.1, ML143345.1, ML143371.1, ML143372.1, ML143377.1, ML143379.1, ML143380.1
    ##   - in 'y': KI270709.1, KI270720.1, KI270721.1, KI270830.1, KI270832.1, KI270869.1, KN196484.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000224.1, GL000251.2, GL000255.2, GL383522.1, GL383556.1, GL383557.1, GL383578.2, GL877875.1, GL949752.1, KI270467.1, KI270706.1, KI270712.1, KI270714.1, KI270735.1, KI270745.1, KI270751.1, KI270765.1, KI270782.1, KI270804.1, KI270809.1, KI270845.1, KI270849.1, KI270850.1, KI270856.1, KI270861.1, KI270868.1, KI270878.1, KI270879.1, KI270897.1, KI270908.1, KI270934.1, KN196479.1, KV575244.1, KV766198.1, KV880764.1, KZ208912.1, KZ208915.1, ML143344.1, ML143345.1, ML143371.1, ML143372.1, ML143377.1, ML143379.1, ML143380.1
    ##   - in 'y': KI270709.1, KI270720.1, KI270721.1, KI270830.1, KI270832.1, KI270869.1, KN196484.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270787.1, KI270803.1, KI270832.1, KV766198.1
    ##   - in 'y': GL000214.1, GL000219.1, GL000220.1, GL000225.1, KI270438.1, KI270442.1, KI270709.1, KI270714.1, KI270728.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270856.1, KI270878.1, KI270879.1, KQ458385.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL339449.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270787.1, KI270803.1, KI270832.1, KV766198.1
    ##   - in 'y': GL000214.1, GL000219.1, GL000220.1, GL000225.1, KI270438.1, KI270442.1, KI270709.1, KI270714.1, KI270728.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270856.1, KI270878.1, KI270879.1, KQ458385.1, ML143372.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270721.1, KI270754.1, KI270819.1, KI270848.1, KI270868.1, KN538370.1, ML143350.1, ML143353.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000219.1, GL000224.1, GL000254.2, GL383522.1, GL383545.1, GL383578.2, GL877875.1, JH159146.1, KI270466.1, KI270467.1, KI270707.1, KI270719.1, KI270731.1, KI270782.1, KI270784.1, KI270850.1, KI270861.1, KI270878.1, KI270879.1, KI270880.1, KI270934.1, KI270936.1, KI270937.1, KV766193.1, KV880768.1, KZ559111.1, ML143366.1, ML143372.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270721.1, KI270754.1, KI270819.1, KI270848.1, KI270868.1, KN538370.1, ML143350.1, ML143353.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000219.1, GL000224.1, GL000254.2, GL383522.1, GL383545.1, GL383578.2, GL877875.1, JH159146.1, KI270466.1, KI270467.1, KI270707.1, KI270719.1, KI270731.1, KI270782.1, KI270784.1, KI270850.1, KI270861.1, KI270878.1, KI270879.1, KI270880.1, KI270934.1, KI270936.1, KI270937.1, KV766193.1, KV880768.1, KZ559111.1, ML143366.1, ML143372.1, ML143375.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KZ208915.1, ML143345.1, ML143355.1
    ##   - in 'y': GL000221.1, KI270711.1, KI270728.1, KI270905.1, KV766198.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KZ208915.1, ML143345.1, ML143355.1
    ##   - in 'y': GL000221.1, KI270711.1, KI270728.1, KI270905.1, KV766198.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL383578.2, KI270538.1, KI270731.1, KI270738.1, KI270751.1, KI270810.1, KI270848.1, KI270861.1, KQ983257.1, KV575244.1, KZ208912.1, KZ559105.1, ML143345.1, ML143367.1, ML143377.1
    ##   - in 'y': GL000218.1, GL000225.1, KI270310.1, KI270337.1, KI270466.1, KI270467.1, KI270721.1, KI270722.1, KI270745.1, KI270821.1, KI270851.1, KI270856.1, KI270880.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL383578.2, KI270538.1, KI270731.1, KI270738.1, KI270751.1, KI270810.1, KI270848.1, KI270861.1, KQ983257.1, KV575244.1, KZ208912.1, KZ559105.1, ML143345.1, ML143367.1, ML143377.1
    ##   - in 'y': GL000218.1, GL000225.1, KI270310.1, KI270337.1, KI270466.1, KI270467.1, KI270721.1, KI270722.1, KI270745.1, KI270821.1, KI270851.1, KI270856.1, KI270880.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270538.1, KI270709.1, KN538372.1, KQ983257.1, ML143345.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000219.1, GL000255.2, GL000256.2, GL339449.2, KI270713.1, KI270714.1, KI270728.1, KI270742.1, KI270816.1, KI270832.1, KI270856.1, KI270860.1, KI270879.1, KI270908.1, KI270937.1, KQ090026.1, KV766196.1, KV766198.1, KV880768.1, ML143367.1, ML143372.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270538.1, KI270709.1, KN538372.1, KQ983257.1, ML143345.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000219.1, GL000255.2, GL000256.2, GL339449.2, KI270713.1, KI270714.1, KI270728.1, KI270742.1, KI270816.1, KI270832.1, KI270856.1, KI270860.1, KI270879.1, KI270908.1, KI270937.1, KQ090026.1, KV766196.1, KV766198.1, KV880768.1, ML143367.1, ML143372.1, ML143377.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383557.1, GL383579.2, KI270784.1, KI270865.1, KN538360.1, KQ458384.1, KV880763.1, ML143352.1, ML143359.1, ML143379.1
    ##   - in 'y': GL000221.1, GL383563.3, GL383581.2, KI270707.1, KI270762.1, KI270809.1, KI270821.1, KI270830.1, KI270844.1, KI270860.1, KI270872.1, KV575243.1, KV575245.1, KZ208914.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383557.1, GL383579.2, KI270784.1, KI270865.1, KN538360.1, KQ458384.1, KV880763.1, ML143352.1, ML143359.1, ML143379.1
    ##   - in 'y': GL000221.1, GL383563.3, GL383581.2, KI270707.1, KI270762.1, KI270809.1, KI270821.1, KI270830.1, KI270844.1, KI270860.1, KI270872.1, KV575243.1, KV575245.1, KZ208914.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000254.2, GL383519.1, GL383578.2, KI270333.1, KI270435.1, KI270442.1, KI270516.1, KI270580.1, KI270583.1, KI270587.1, KI270706.1, KI270718.1, KI270719.1, KI270725.1, KI270730.1, KI270736.1, KI270749.1, KI270754.1, KI270765.1, KI270772.1, KI270831.1, KI270849.1, KQ031384.1, ML143352.1, ML143371.1
    ##   - in 'y': GL000214.1, GL000221.1, GL000251.2, KI270709.1, KI270721.1, KI270724.1, KI270742.1, KI270757.1, KI270830.1, KI270844.1, KI270857.1, KQ031389.1, KV880768.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000254.2, GL383519.1, GL383578.2, KI270333.1, KI270435.1, KI270442.1, KI270516.1, KI270580.1, KI270583.1, KI270587.1, KI270706.1, KI270718.1, KI270719.1, KI270725.1, KI270730.1, KI270736.1, KI270749.1, KI270754.1, KI270765.1, KI270772.1, KI270831.1, KI270849.1, KQ031384.1, ML143352.1, ML143371.1
    ##   - in 'y': GL000214.1, GL000221.1, GL000251.2, KI270709.1, KI270721.1, KI270724.1, KI270742.1, KI270757.1, KI270830.1, KI270844.1, KI270857.1, KQ031389.1, KV880768.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000251.2, GL000254.2, GL949746.1, GL949752.1, KI270442.1, KI270707.1, KI270822.1, KI270827.1, KI270848.1, KI270850.1, KI270868.1, KI270872.1, KI270876.1, KI270896.1, KN196477.1, KN538361.1, KN538364.1, KN538370.1, KQ031389.1, KV575244.1, KV880764.1, KZ208907.1, ML143353.1, ML143355.1, ML143371.1, ML143375.1, ML143377.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, KI270467.1, KI270856.1, KQ458385.1, ML143382.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000251.2, GL000254.2, GL949746.1, GL949752.1, KI270442.1, KI270707.1, KI270822.1, KI270827.1, KI270848.1, KI270850.1, KI270868.1, KI270872.1, KI270876.1, KI270896.1, KN196477.1, KN538361.1, KN538364.1, KN538370.1, KQ031389.1, KV575244.1, KV880764.1, KZ208907.1, ML143353.1, ML143355.1, ML143371.1, ML143375.1, ML143377.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, KI270467.1, KI270856.1, KQ458385.1, ML143382.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270330.1, KI270538.1, KI270714.1, KI270720.1, KI270738.1, KI270765.1, KN538370.1, KZ208915.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000256.2, KI270310.1, KI270337.1, KI270466.1, KI270467.1, KI270742.1, KI270744.1, KI270754.1, KQ090026.1, ML143345.1, ML143366.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000224.1, KI270330.1, KI270538.1, KI270714.1, KI270720.1, KI270738.1, KI270765.1, KN538370.1, KZ208915.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000256.2, KI270310.1, KI270337.1, KI270466.1, KI270467.1, KI270742.1, KI270744.1, KI270754.1, KQ090026.1, ML143345.1, ML143366.1, ML143372.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270728.1, KI270745.1, KI270784.1, KI270848.1, KI270850.1, KQ458383.1, KV766198.1, KZ208907.1, ML143344.1, ML143352.1
    ##   - in 'y': GL000224.1, KI270330.1, KI270438.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270728.1, KI270745.1, KI270784.1, KI270848.1, KI270850.1, KQ458383.1, KV766198.1, KZ208907.1, ML143344.1, ML143352.1
    ##   - in 'y': GL000224.1, KI270330.1, KI270438.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270711.1, KI270742.1, KI270745.1, KI270784.1, KI270803.1, KI270830.1, KI270848.1, KI270849.1, KI270850.1, KI270851.1, KI270857.1, KQ458383.1, KV766198.1, KZ208907.1, KZ208912.1, ML143345.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL383563.3, GL383577.2, KI270320.1, KI270435.1, KI270508.1, KI270517.1, KI270519.1, KI270538.1, KI270587.1, KI270589.1, KI270591.1, KI270712.1, KI270718.1, KI270724.1, KI270725.1, KI270729.1, KI270730.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270746.1, KI270750.1, KI270751.1, KI270754.1, KI270756.1, KI270757.1, KI270772.1, KI270805.1, KI270861.1, KI270869.1, KI270936.1, KN196487.1, KN538364.1, KN538372.1, KQ031384.1, KQ031389.1, KQ983257.1, KV880764.1, KV880768.1, KZ208916.1, ML143341.1, ML143353.1, ML143354.1, ML143362.1, ML143364.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270711.1, KI270742.1, KI270745.1, KI270784.1, KI270803.1, KI270830.1, KI270848.1, KI270849.1, KI270850.1, KI270851.1, KI270857.1, KQ458383.1, KV766198.1, KZ208907.1, KZ208912.1, ML143345.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000008.2, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL383563.3, GL383577.2, KI270320.1, KI270435.1, KI270508.1, KI270517.1, KI270519.1, KI270538.1, KI270587.1, KI270589.1, KI270591.1, KI270712.1, KI270718.1, KI270724.1, KI270725.1, KI270729.1, KI270730.1, KI270733.1, KI270734.1, KI270735.1, KI270736.1, KI270746.1, KI270750.1, KI270751.1, KI270754.1, KI270756.1, KI270757.1, KI270772.1, KI270805.1, KI270861.1, KI270869.1, KI270936.1, KN196487.1, KN538364.1, KN538372.1, KQ031384.1, KQ031389.1, KQ983257.1, KV880764.1, KV880768.1, KZ208916.1, ML143341.1, ML143353.1, ML143354.1, ML143362.1, ML143364.1, ML143365.1, ML143370.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270711.1, KI270742.1, KI270783.1, KI270784.1, KI270792.1, KI270803.1, KI270813.1, KI270848.1, KI270849.1, KI270850.1, KI270853.1, KQ458383.1, KV766198.1, KZ208907.1, KZ208912.1, ML143344.1, ML143366.1, ML143371.1, ML143377.1, KI270709.1, GL000218.1, GL000220.1, GL000225.1, GL383563.3, GL383577.2, KI270320.1, KI270435.1, KI270508.1, KI270517.1, KI270519.1, KI270538.1, KI270587.1, KI270589.1, KI270591.1, KI270724.1, KI270725.1, KI270734.1, KI270746.1, KI270750.1, KI270751.1, KI270756.1, KI270772.1, KI270805.1, KI270861.1, KI270936.1, KN196487.1, KN538372.1, KQ983257.1, KV880768.1, KZ208916.1, ML143341.1, ML143354.1, ML143364.1, ML143370.1
    ##   - in 'y': ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, GL000251.2, GL000252.2, GL000254.2, GL383578.2, KI270711.1, KI270742.1, KI270783.1, KI270784.1, KI270792.1, KI270803.1, KI270813.1, KI270848.1, KI270849.1, KI270850.1, KI270853.1, KQ458383.1, KV766198.1, KZ208907.1, KZ208912.1, ML143344.1, ML143366.1, ML143371.1, ML143377.1, KI270709.1, GL000218.1, GL000220.1, GL000225.1, GL383563.3, GL383577.2, KI270320.1, KI270435.1, KI270508.1, KI270517.1, KI270519.1, KI270538.1, KI270587.1, KI270589.1, KI270591.1, KI270724.1, KI270725.1, KI270734.1, KI270746.1, KI270750.1, KI270751.1, KI270756.1, KI270772.1, KI270805.1, KI270861.1, KI270936.1, KN196487.1, KN538372.1, KQ983257.1, KV880768.1, KZ208916.1, ML143341.1, ML143354.1, ML143364.1, ML143370.1
    ##   - in 'y': ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, chrY
    ##   - in 'y': GL000224.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, chrY
    ##   - in 'y': GL000224.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270853.1, chrY, GL000224.1
    ##   - in 'y': GL000194.1, GL000251.2, GL000256.2, KI270438.1, KI270713.1, KI270728.1, KI270742.1, KI270744.1, KI270765.1, KI270782.1, KI270871.1, KI270899.1, KV575244.1, KV766198.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270853.1, chrY, GL000224.1
    ##   - in 'y': GL000194.1, GL000251.2, GL000256.2, KI270438.1, KI270713.1, KI270728.1, KI270742.1, KI270744.1, KI270765.1, KI270782.1, KI270871.1, KI270899.1, KV575244.1, KV766198.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, ML143380.1, chrY, KI270728.1, KI270765.1, KI270871.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000225.1, GL000255.2, GL383522.1, GL949752.1, KI270711.1, KI270714.1, KI270733.1, KI270734.1, KI270736.1, KI270745.1, KI270754.1, KI270784.1, KI270880.1, KI270897.1, KI270908.1, KN196487.1, KQ458383.1, KZ208915.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, ML143380.1, chrY, KI270728.1, KI270765.1, KI270871.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000225.1, GL000255.2, GL383522.1, GL949752.1, KI270711.1, KI270714.1, KI270733.1, KI270734.1, KI270736.1, KI270745.1, KI270754.1, KI270784.1, KI270880.1, KI270897.1, KI270908.1, KN196487.1, KQ458383.1, KZ208915.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270765.1, KI270871.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL383522.1, GL949752.1, KI270711.1, KI270714.1, KI270733.1, KI270734.1, KI270736.1, KI270745.1, KI270754.1, KI270784.1, KI270853.1, KI270880.1, KI270897.1, KI270908.1, KN196487.1, KQ458383.1, KZ208915.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270765.1, KI270871.1, KI270899.1, KV880764.1, ML143352.1, ML143355.1, ML143359.1, ML143375.1, ML143379.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL383522.1, GL949752.1, KI270711.1, KI270714.1, KI270733.1, KI270734.1, KI270736.1, KI270745.1, KI270754.1, KI270784.1, KI270853.1, KI270880.1, KI270897.1, KI270908.1, KN196487.1, KQ458383.1, KZ208915.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270712.1, KI270745.1, KI270782.1, KI270803.1, KN538370.1, ML143345.1, ML143353.1, ML143355.1, ML143358.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000219.1, GL000221.1, GL000225.1, GL000253.2, GL339449.2, GL383563.3, GL949752.1, JH159146.1, KI270442.1, KI270706.1, KI270711.1, KI270721.1, KI270731.1, KI270734.1, KI270751.1, KI270783.1, KI270816.1, KI270817.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270844.1, KI270856.1, KI270861.1, KI270897.1, KI270905.1, KN196484.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KZ208913.1, KZ559105.1, ML143367.1, ML143371.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270712.1, KI270745.1, KI270782.1, KI270803.1, KN538370.1, ML143345.1, ML143353.1, ML143355.1, ML143358.1, ML143375.1, ML143379.1
    ##   - in 'y': GL000008.2, GL000214.1, GL000216.2, GL000219.1, GL000221.1, GL000225.1, GL000253.2, GL339449.2, GL383563.3, GL949752.1, JH159146.1, KI270442.1, KI270706.1, KI270711.1, KI270721.1, KI270731.1, KI270734.1, KI270751.1, KI270783.1, KI270816.1, KI270817.1, KI270819.1, KI270821.1, KI270830.1, KI270831.1, KI270844.1, KI270856.1, KI270861.1, KI270897.1, KI270905.1, KN196484.1, KN538364.1, KN538372.1, KQ031389.1, KQ090026.1, KZ208913.1, KZ559105.1, ML143367.1, ML143371.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, GL383575.2, JH159146.1, KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270729.1, KI270732.1, KI270733.1, KI270821.1, KI270876.1, KI270894.1, KN196478.1, KQ031384.1, KQ031389.1, KV880764.1, KZ208913.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000218.1, GL000220.1, GL000251.2, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270832.1, KI270844.1, KI270851.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KI270899.1, KN196484.1, KQ458383.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, KZ559112.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383563.3, GL383575.2, JH159146.1, KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270729.1, KI270732.1, KI270733.1, KI270821.1, KI270876.1, KI270894.1, KN196478.1, KQ031384.1, KQ031389.1, KV880764.1, KZ208913.1, ML143352.1, ML143355.1
    ##   - in 'y': GL000009.2, GL000218.1, GL000220.1, GL000251.2, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270832.1, KI270844.1, KI270851.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KI270899.1, KN196484.1, KQ458383.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, KZ559112.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000253.2, GL383563.3, GL383575.2, GL383578.2, JH159146.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270714.1, KI270729.1, KI270742.1, KI270821.1, KI270831.1, KI270850.1, KI270853.1, KI270876.1, KI270879.1, KI270908.1, KN196478.1, KQ031384.1, KQ031389.1, KZ208913.1, KZ559103.1, ML143345.1, ML143352.1, ML143371.1, ML143380.1, GL000009.2, GL000218.1, GL000251.2, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270851.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KN196484.1, KQ458383.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, KZ559112.1, ML143377.1
    ##   - in 'y': KI270330.1, KI270731.1, KI270869.1, KN196487.1, KN538370.1, KZ208912.1, ML143353.1, ML143358.1, ML143375.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000253.2, GL383563.3, GL383575.2, GL383578.2, JH159146.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270713.1, KI270714.1, KI270729.1, KI270742.1, KI270821.1, KI270831.1, KI270850.1, KI270853.1, KI270876.1, KI270879.1, KI270908.1, KN196478.1, KQ031384.1, KQ031389.1, KZ208913.1, KZ559103.1, ML143345.1, ML143352.1, ML143371.1, ML143380.1, GL000009.2, GL000218.1, GL000251.2, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270782.1, KI270783.1, KI270816.1, KI270830.1, KI270851.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KN196484.1, KQ458383.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, KZ559112.1, ML143377.1
    ##   - in 'y': KI270330.1, KI270731.1, KI270869.1, KN196487.1, KN538370.1, KZ208912.1, ML143353.1, ML143358.1, ML143375.1, ML143378.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL339449.2, GL383563.3, GL383575.2, GL383578.2, GL877875.1, JH159146.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270803.1, KI270821.1, KI270849.1, KI270876.1, KI270894.1, KI270908.1, KN196478.1, KQ031384.1, KQ031389.1, KZ208913.1, KZ559103.1, GL000009.2, GL000218.1, GL000220.1, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270783.1, KI270816.1, KI270830.1, KI270844.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KN196484.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, ML143377.1, KI270330.1, KI270731.1, KI270869.1, KN196487.1, ML143358.1, ML143378.1
    ##   - in 'y': KI270718.1, KN196479.1, KQ090016.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000219.1, GL000224.1, GL000225.1, GL339449.2, GL383563.3, GL383575.2, GL383578.2, GL877875.1, JH159146.1, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270729.1, KI270733.1, KI270803.1, KI270821.1, KI270849.1, KI270876.1, KI270894.1, KI270908.1, KN196478.1, KQ031384.1, KQ031389.1, KZ208913.1, KZ559103.1, GL000009.2, GL000218.1, GL000220.1, GL383522.1, GL383574.1, KI270582.1, KI270743.1, KI270762.1, KI270783.1, KI270816.1, KI270830.1, KI270844.1, KI270856.1, KI270857.1, KI270865.1, KI270868.1, KI270878.1, KI270897.1, KN196484.1, KQ983257.1, KV575244.1, KV880763.1, KZ559111.1, ML143377.1, KI270330.1, KI270731.1, KI270869.1, KN196487.1, ML143358.1, ML143378.1
    ##   - in 'y': KI270718.1, KN196479.1, KQ090016.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, ML143380.1, chr18
    ##   - in 'y': GL000216.2, GL000220.1, KI270713.1, chr13, chr20, chr21
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270712.1, ML143380.1, chr18
    ##   - in 'y': GL000216.2, GL000220.1, KI270713.1, chr13, chr20, chr21
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, ML143380.1, GL000216.2, GL000220.1
    ##   - in 'y': GL000253.2, GL000256.2, KI270711.1, KI270714.1, KI270742.1, KI270765.1, KI270844.1, KI270853.1, KI270857.1, KI270861.1, KQ090026.1, KV575244.1, KV766198.1, KZ559112.1, ML143345.1, ML143352.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, ML143380.1, GL000216.2, GL000220.1
    ##   - in 'y': GL000253.2, GL000256.2, KI270711.1, KI270714.1, KI270742.1, KI270765.1, KI270844.1, KI270853.1, KI270857.1, KI270861.1, KQ090026.1, KV575244.1, KV766198.1, KZ559112.1, ML143345.1, ML143352.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383575.2, GL383578.2, GL949747.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270731.1, KI270745.1, KI270765.1, KI270831.1, KI270850.1, KI270853.1, KQ031389.1, KV766198.1, ML143371.1, ML143372.1
    ##   - in 'y': GL000009.2, GL000214.1, GL877875.1, KI270438.1, KI270538.1, KI270707.1, KI270708.1, KI270733.1, KI270816.1, KI270844.1, KI270905.1, KN196484.1, KV575244.1, KV880768.1, KZ559103.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000251.2, GL383575.2, GL383578.2, GL949747.2, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270731.1, KI270745.1, KI270765.1, KI270831.1, KI270850.1, KI270853.1, KQ031389.1, KV766198.1, ML143371.1, ML143372.1
    ##   - in 'y': GL000009.2, GL000214.1, GL877875.1, KI270438.1, KI270538.1, KI270707.1, KI270708.1, KI270733.1, KI270816.1, KI270844.1, KI270905.1, KN196484.1, KV575244.1, KV880768.1, KZ559103.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270438.1, KI270509.1, KI270709.1, KI270742.1, KI270744.1, KI270880.1, KN538364.1, KQ983257.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': chrM
    ##   - in 'y': GL000214.1, GL000216.2, GL000219.1, GL000220.1, GL000224.1, GL000225.1, KI270330.1, KI270438.1, KI270509.1, KI270709.1, KI270742.1, KI270744.1, KI270880.1, KN538364.1, KQ983257.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270467.1
    ##   - in 'y': GL000214.1, GL000225.1, KI270310.1, KI270538.1, KI270731.1, KN196487.1, KQ031384.1, ML143372.1, ML143380.1, chr18, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270467.1
    ##   - in 'y': GL000214.1, GL000225.1, KI270310.1, KI270538.1, KI270731.1, KN196487.1, KQ031384.1, ML143372.1, ML143380.1, chr18, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, KI270438.1, KI270709.1, KI270765.1, KQ031389.1, KV880764.1, ML143355.1, ML143379.1, chr21, chr5
    ##   - in 'y': KI270712.1, chr14, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, KI270438.1, KI270709.1, KI270765.1, KQ031389.1, KV880764.1, ML143355.1, ML143379.1, chr21, chr5
    ##   - in 'y': KI270712.1, chr14, chrX
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': ML143352.1
    ##   - in 'y': GL383563.3, KI270330.1, chr18
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': ML143352.1
    ##   - in 'y': GL383563.3, KI270330.1, chr18
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949746.1, KI270709.1, KI270712.1, KI270728.1, KI270765.1, KI270899.1, KV880764.1, ML143355.1, ML143375.1
    ##   - in 'y': GL000221.1, GL000255.2, KI270442.1, KI270714.1, KI270744.1, KI270783.1, KI270807.1, KI270816.1, KI270819.1, KI270830.1, KI270831.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270905.1, KN196484.1, KQ090026.1, KQ458383.1, KV575244.1, KV766198.1, KV880768.1, KZ208908.1, KZ559112.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL949746.1, KI270709.1, KI270712.1, KI270728.1, KI270765.1, KI270899.1, KV880764.1, ML143355.1, ML143375.1
    ##   - in 'y': GL000221.1, GL000255.2, KI270442.1, KI270714.1, KI270744.1, KI270783.1, KI270807.1, KI270816.1, KI270819.1, KI270830.1, KI270831.1, KI270849.1, KI270850.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270879.1, KI270880.1, KI270897.1, KI270905.1, KN196484.1, KQ090026.1, KQ458383.1, KV575244.1, KV766198.1, KV880768.1, KZ208908.1, KZ559112.1, ML143371.1, ML143372.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383520.2, GL383575.2, KI270733.1, KI270859.1, KN538372.1, KQ458385.1, KV766196.1
    ##   - in 'y': GL000218.1, GL383578.2, GL383579.2, GL383580.2, GL949752.1, KI270582.1, KI270709.1, KI270712.1, KI270719.1, KI270732.1, KI270749.1, KI270762.1, KI270772.1, KI270792.1, KI270824.1, KI270840.1, KI270849.1, KI270855.1, KI270861.1, KI270870.1, KI270871.1, KI270894.1, KI270907.1, KI270925.1, KN196484.1, KN538370.1, KQ031389.1, KQ090026.1, KQ458383.1, KQ458384.1, KV575243.1, KZ208906.1, KZ208912.1, KZ559109.1, ML143352.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL383520.2, GL383575.2, KI270733.1, KI270859.1, KN538372.1, KQ458385.1, KV766196.1
    ##   - in 'y': GL000218.1, GL383578.2, GL383579.2, GL383580.2, GL949752.1, KI270582.1, KI270709.1, KI270712.1, KI270719.1, KI270732.1, KI270749.1, KI270762.1, KI270772.1, KI270792.1, KI270824.1, KI270840.1, KI270849.1, KI270855.1, KI270861.1, KI270870.1, KI270871.1, KI270894.1, KI270907.1, KI270925.1, KN196484.1, KN538370.1, KQ031389.1, KQ090026.1, KQ458383.1, KQ458384.1, KV575243.1, KZ208906.1, KZ208912.1, KZ559109.1, ML143352.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL383519.1, KI270734.1, KI270750.1, KI270754.1, KI270765.1, KV880764.1, KZ208922.1, KZ559105.1, ML143352.1, ML143355.1, ML143369.1
    ##   - in 'y': GL000214.1, GL000224.1, GL383526.1, GL383556.1, GL383578.2, KI270707.1, KI270731.1, KI270770.1, KI270787.1, KI270809.1, KI270819.1, KI270821.1, KI270831.1, KI270836.1, KI270904.1, KN196484.1, KQ090026.1, KQ458383.1, KV880763.1, KZ208912.1, ML143353.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000205.2, GL383519.1, KI270734.1, KI270750.1, KI270754.1, KI270765.1, KV880764.1, KZ208922.1, KZ559105.1, ML143352.1, ML143355.1, ML143369.1
    ##   - in 'y': GL000214.1, GL000224.1, GL383526.1, GL383556.1, GL383578.2, KI270707.1, KI270731.1, KI270770.1, KI270787.1, KI270809.1, KI270819.1, KI270821.1, KI270831.1, KI270836.1, KI270904.1, KN196484.1, KQ090026.1, KQ458383.1, KV880763.1, KZ208912.1, ML143353.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270435.1, KI270589.1, KI270718.1, KI270733.1, KI270836.1, KQ031384.1, KV766198.1, ML143344.1, ML143371.1
    ##   - in 'y': GL000219.1, GL000225.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270435.1, KI270589.1, KI270718.1, KI270733.1, KI270836.1, KQ031384.1, KV766198.1, ML143344.1, ML143371.1
    ##   - in 'y': GL000219.1, GL000225.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000216.2, GL000220.1, GL000255.2, KI270435.1, KI270589.1, KI270718.1, KI270730.1, KI270736.1, KI270742.1, KI270751.1, KI270757.1, KI270836.1, KI270902.1, KI270924.1, KQ031384.1, KV766198.1, ML143344.1, ML143345.1, ML143362.1, ML143370.1, ML143371.1, GL000219.1, GL000225.1, KI270712.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1
    ##   - in 'y': GL000194.1, KI270467.1, KI270731.1, KN538364.1, KQ031389.1, KZ208915.1, KZ559112.1, ML143353.1, ML143355.1, ML143365.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000216.2, GL000220.1, GL000255.2, KI270435.1, KI270589.1, KI270718.1, KI270730.1, KI270736.1, KI270742.1, KI270751.1, KI270757.1, KI270836.1, KI270902.1, KI270924.1, KQ031384.1, KV766198.1, ML143344.1, ML143345.1, ML143362.1, ML143370.1, ML143371.1, GL000219.1, GL000225.1, KI270712.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1
    ##   - in 'y': GL000194.1, KI270467.1, KI270731.1, KN538364.1, KQ031389.1, KZ208915.1, KZ559112.1, ML143353.1, ML143355.1, ML143365.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270337.1, KI270709.1, KI270742.1, KI270757.1, KI270836.1, KI270902.1, KI270924.1, KV766198.1, ML143344.1, ML143345.1, ML143362.1, ML143370.1, GL000219.1, GL000225.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1, GL000194.1, KI270467.1, KI270731.1, KN538364.1, ML143353.1, ML143375.1
    ##   - in 'y': GL000253.2, KI270519.1, KI270735.1, KI270749.1, KI270754.1, KI270765.1, KV880764.1, ML143352.1, ML143360.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270337.1, KI270709.1, KI270742.1, KI270757.1, KI270836.1, KI270902.1, KI270924.1, KV766198.1, ML143344.1, ML143345.1, ML143362.1, ML143370.1, GL000219.1, GL000225.1, KI270713.1, KI270714.1, KI270738.1, KI270762.1, KI270779.1, KI270894.1, KQ983257.1, KZ208912.1, ML143377.1, GL000194.1, KI270467.1, KI270731.1, KN538364.1, ML143353.1, ML143375.1
    ##   - in 'y': GL000253.2, KI270519.1, KI270735.1, KI270749.1, KI270754.1, KI270765.1, KV880764.1, ML143352.1, ML143360.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270742.1, KI270765.1, ML143380.1, chr13, chr18, chr22
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270442.1, KI270729.1, KI270733.1, KI270744.1, KN196487.1, chr9
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270742.1, KI270765.1, ML143380.1, chr13, chr18, chr22
    ##   - in 'y': GL000214.1, GL000216.2, GL000224.1, GL000225.1, KI270442.1, KI270729.1, KI270733.1, KI270744.1, KN196487.1, chr9
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270742.1, ML143380.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000225.1, GL000255.2, KI270330.1, KI270709.1, KI270728.1, KI270733.1, KI270736.1, KI270744.1, KI270816.1, KI270880.1, KI270937.1, KN196484.1, KN196487.1, KZ208915.1, ML143345.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270742.1, ML143380.1
    ##   - in 'y': GL000194.1, GL000214.1, GL000216.2, GL000219.1, GL000225.1, GL000255.2, KI270330.1, KI270709.1, KI270728.1, KI270733.1, KI270736.1, KI270744.1, KI270816.1, KI270880.1, KI270937.1, KN196484.1, KN196487.1, KZ208915.1, ML143345.1, ML143377.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270438.1, GL000194.1, GL000216.2, GL000225.1, KI270330.1, KI270709.1, KI270733.1, KI270736.1, KI270816.1, KI270937.1, KN196484.1, KN196487.1
    ##   - in 'y': GL000251.2, GL000253.2, GL339449.2, GL949746.1, KI270707.1, KI270711.1, KI270718.1, KI270765.1, KI270784.1, KI270830.1, KI270849.1, KI270850.1, KI270856.1, KI270897.1, KI270905.1, KI270908.1, KN538364.1, KN538370.1, KQ031389.1, KQ090026.1, KQ458383.1, KV575244.1, KV880764.1, KZ559112.1, ML143352.1, ML143353.1, ML143355.1, ML143371.1, ML143372.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, KI270438.1, GL000194.1, GL000216.2, GL000225.1, KI270330.1, KI270709.1, KI270733.1, KI270736.1, KI270816.1, KI270937.1, KN196484.1, KN196487.1
    ##   - in 'y': GL000251.2, GL000253.2, GL339449.2, GL949746.1, KI270707.1, KI270711.1, KI270718.1, KI270765.1, KI270784.1, KI270830.1, KI270849.1, KI270850.1, KI270856.1, KI270897.1, KI270905.1, KI270908.1, KN538364.1, KN538370.1, KQ031389.1, KQ090026.1, KQ458383.1, KV575244.1, KV880764.1, KZ559112.1, ML143352.1, ML143353.1, ML143355.1, ML143371.1, ML143372.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, JH159146.1, KI270712.1, KI270734.1, KI270738.1, KI270779.1, KI270803.1, KI270819.1, KI270838.1, KI270848.1, KI270853.1, KI270905.1, KI270926.1, KI270937.1, KN538361.1, KV880764.1, ML143355.1, ML143358.1, ML143366.1, ML143367.1, ML143372.1, ML143377.1
    ##   - in 'y': GL000009.2, GL000251.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KV880768.1, KZ559103.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, JH159146.1, KI270712.1, KI270734.1, KI270738.1, KI270779.1, KI270803.1, KI270819.1, KI270838.1, KI270848.1, KI270853.1, KI270905.1, KI270926.1, KI270937.1, KN538361.1, KV880764.1, ML143355.1, ML143358.1, ML143366.1, ML143367.1, ML143372.1, ML143377.1
    ##   - in 'y': GL000009.2, GL000251.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KV880768.1, KZ559103.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, JH159146.1, KI270442.1, KI270712.1, KI270728.1, KI270738.1, KI270802.1, KI270803.1, KI270819.1, KI270838.1, KI270848.1, KI270868.1, KI270926.1, KI270937.1, KN538361.1, ML143352.1, ML143358.1, ML143366.1, GL000009.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KZ559103.1, chrM
    ##   - in 'y': GL383557.1, GL383563.3, GL877875.1, GL949746.1, KI270784.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, JH159146.1, KI270442.1, KI270712.1, KI270728.1, KI270738.1, KI270802.1, KI270803.1, KI270819.1, KI270838.1, KI270848.1, KI270868.1, KI270926.1, KI270937.1, KN538361.1, ML143352.1, ML143358.1, ML143366.1, GL000009.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KZ559103.1, chrM
    ##   - in 'y': GL383557.1, GL383563.3, GL877875.1, GL949746.1, KI270784.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, JH159146.1, KI270712.1, KI270734.1, KI270738.1, KI270802.1, KI270803.1, KI270819.1, KI270838.1, KI270868.1, KI270926.1, KI270937.1, KN538361.1, KQ458385.1, ML143358.1, ML143367.1, ML143377.1, GL000009.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KZ559103.1, GL383563.3, GL877875.1, KI270784.1
    ##   - in 'y': GL383522.1, GL949752.1, KI270745.1, KN538364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, JH159146.1, KI270712.1, KI270734.1, KI270738.1, KI270802.1, KI270803.1, KI270819.1, KI270838.1, KI270868.1, KI270926.1, KI270937.1, KN538361.1, KQ458385.1, ML143358.1, ML143367.1, ML143377.1, GL000009.2, GL383578.2, KI270337.1, KI270467.1, KI270707.1, KI270821.1, KZ559103.1, GL383563.3, GL877875.1, KI270784.1
    ##   - in 'y': GL383522.1, GL949752.1, KI270745.1, KN538364.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270467.1, KI270709.1, KI270712.1, KI270720.1, KI270738.1, KI270784.1, KI270803.1, KV766198.1, KZ559112.1
    ##   - in 'y': GL000009.2, GL000194.1, KI270706.1, KI270725.1, KI270816.1, KI270819.1, KI270830.1, KI270836.1, KI270853.1, KI270897.1, KV575244.1, KZ559105.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270467.1, KI270709.1, KI270712.1, KI270720.1, KI270738.1, KI270784.1, KI270803.1, KV766198.1, KZ559112.1
    ##   - in 'y': GL000009.2, GL000194.1, KI270706.1, KI270725.1, KI270816.1, KI270819.1, KI270830.1, KI270836.1, KI270853.1, KI270897.1, KV575244.1, KZ559105.1, ML143372.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000220.1, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270733.1, KI270765.1, KI270799.1, KI270831.1, KI270892.1, KI270905.1, KI270907.1, KZ208913.1, ML143345.1, ML143355.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000255.2, GL949747.2, KI270709.1, KI270714.1, KI270745.1, KI270761.1, KI270782.1, KI270805.1, KI270808.1, KI270816.1, KI270830.1, KI270836.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KI270880.1, KI270924.1, KN196479.1, KN196484.1, KQ031389.1, KQ759759.1, KV880768.1, ML143353.1, ML143382.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000220.1, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270728.1, KI270733.1, KI270765.1, KI270799.1, KI270831.1, KI270892.1, KI270905.1, KI270907.1, KZ208913.1, ML143345.1, ML143355.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000255.2, GL949747.2, KI270709.1, KI270714.1, KI270745.1, KI270761.1, KI270782.1, KI270805.1, KI270808.1, KI270816.1, KI270830.1, KI270836.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KI270880.1, KI270924.1, KN196479.1, KN196484.1, KQ031389.1, KQ759759.1, KV880768.1, ML143353.1, ML143382.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000219.1, GL383578.2, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270765.1, KI270799.1, KI270822.1, KI270831.1, KI270851.1, KI270853.1, KI270892.1, KI270905.1, KI270907.1, KV575244.1, KZ208913.1, ML143371.1, ML143372.1, ML143380.1, GL000255.2, KI270745.1, KI270761.1, KI270782.1, KI270808.1, KI270816.1, KI270830.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KN196484.1, KQ759759.1, ML143382.1
    ##   - in 'y': KI270768.1, KN538364.1, KV880764.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000219.1, GL383578.2, GL877875.1, KI270337.1, KI270466.1, KI270467.1, KI270706.1, KI270728.1, KI270765.1, KI270799.1, KI270822.1, KI270831.1, KI270851.1, KI270853.1, KI270892.1, KI270905.1, KI270907.1, KV575244.1, KZ208913.1, ML143371.1, ML143372.1, ML143380.1, GL000255.2, KI270745.1, KI270761.1, KI270782.1, KI270808.1, KI270816.1, KI270830.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KN196484.1, KQ759759.1, ML143382.1
    ##   - in 'y': KI270768.1, KN538364.1, KV880764.1, ML143375.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000219.1, GL000253.2, GL877875.1, KI270466.1, KI270467.1, KI270799.1, KI270822.1, KI270831.1, KI270851.1, KI270853.1, KI270861.1, KI270892.1, KI270905.1, KI270907.1, KV575244.1, KZ208913.1, GL949747.2, KI270709.1, KI270714.1, KI270745.1, KI270761.1, KI270805.1, KI270808.1, KI270816.1, KI270830.1, KI270836.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KI270924.1, KN196479.1, KN196484.1, KQ759759.1, KV880768.1, ML143382.1, KI270768.1
    ##   - in 'y': GL000224.1, GL000251.2, GL949746.1, KI270330.1, KI270435.1, KI270711.1, KI270718.1, KI270736.1, KI270754.1, KI270857.1, KI270899.1, KN538370.1, KQ090023.1, ML143350.1, ML143358.1, ML143359.1, ML143362.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000009.2, GL000219.1, GL000253.2, GL877875.1, KI270466.1, KI270467.1, KI270799.1, KI270822.1, KI270831.1, KI270851.1, KI270853.1, KI270861.1, KI270892.1, KI270905.1, KI270907.1, KV575244.1, KZ208913.1, GL949747.2, KI270709.1, KI270714.1, KI270745.1, KI270761.1, KI270805.1, KI270808.1, KI270816.1, KI270830.1, KI270836.1, KI270844.1, KI270856.1, KI270876.1, KI270879.1, KI270924.1, KN196479.1, KN196484.1, KQ759759.1, KV880768.1, ML143382.1, KI270768.1
    ##   - in 'y': GL000224.1, GL000251.2, GL949746.1, KI270330.1, KI270435.1, KI270711.1, KI270718.1, KI270736.1, KI270754.1, KI270857.1, KI270899.1, KN538370.1, KQ090023.1, ML143350.1, ML143358.1, ML143359.1, ML143362.1, ML143365.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270803.1, KI270819.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1
    ##   - in 'y': GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270803.1, KI270819.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270714.1, KI270733.1, KI270744.1, KI270856.1, ML143352.1, ML143365.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270783.1, KI270849.1, KZ559112.1
    ##   - in 'y': GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KQ031389.1, KV880764.1, KZ559109.1, ML143344.1, ML143345.1, ML143350.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270714.1, KI270733.1, KI270744.1, KI270856.1, ML143352.1, ML143365.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270783.1, KI270849.1, KZ559112.1
    ##   - in 'y': GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KQ031389.1, KV880764.1, KZ559109.1, ML143344.1, ML143345.1, ML143350.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270733.1, KI270744.1, KI270745.1, KI270818.1, KI270856.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270849.1, KI270857.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1, GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KV880764.1, KZ559109.1, ML143344.1, ML143345.1, ML143350.1, ML143379.1
    ##   - in 'y': GL000216.2, KI270879.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270733.1, KI270744.1, KI270745.1, KI270818.1, KI270856.1, ML143352.1, ML143353.1, ML143355.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, GL949752.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270802.1, KI270849.1, KI270857.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1, GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KV880764.1, KZ559109.1, ML143344.1, ML143345.1, ML143350.1, ML143379.1
    ##   - in 'y': GL000216.2, KI270879.1, ML143366.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, KI270310.1, KI270714.1, KI270745.1, ML143365.1, ML143372.1, GL000009.2, GL000218.1, GL949752.1, KI270754.1, KI270782.1, KN538361.1, GL383526.1, KI270784.1, KI270895.1, KI270897.1, KI270903.1, KZ559109.1, ML143350.1, KI270879.1, ML143366.1
    ##   - in 'y': GL000008.2, GL000224.1, GL000254.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270830.1, KI270848.1, KI270850.1, KI270862.1, KI270866.1, KQ031384.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, KI270310.1, KI270714.1, KI270745.1, ML143365.1, ML143372.1, GL000009.2, GL000218.1, GL949752.1, KI270754.1, KI270782.1, KN538361.1, GL383526.1, KI270784.1, KI270895.1, KI270897.1, KI270903.1, KZ559109.1, ML143350.1, KI270879.1, ML143366.1
    ##   - in 'y': GL000008.2, GL000224.1, GL000254.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270830.1, KI270848.1, KI270850.1, KI270862.1, KI270866.1, KQ031384.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270706.1, KI270744.1, KI270745.1, KI270818.1, KI270905.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270803.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1, GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KQ031389.1, KZ559109.1, ML143344.1, ML143350.1, ML143379.1, GL000216.2, KI270879.1, ML143366.1, GL000008.2, GL383520.2, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270830.1, KI270848.1, KI270850.1, KI270862.1, KI270866.1, KQ031384.1, KV575244.1
    ##   - in 'y': GL383578.2, KI270707.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000225.1, GL000251.2, KI270310.1, KI270706.1, KI270744.1, KI270745.1, KI270818.1, KI270905.1, KN538364.1, KV766198.1, ML143352.1, ML143353.1, ML143355.1, ML143364.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270754.1, KI270782.1, KI270783.1, KI270803.1, KI270849.1, KI270853.1, KI270857.1, KI270878.1, KI270894.1, KN538361.1, KZ559112.1, ML143375.1, ML143377.1, ML143380.1, GL000253.2, GL383526.1, GL949746.1, KI270784.1, KI270831.1, KI270868.1, KI270880.1, KI270895.1, KI270897.1, KI270903.1, KN538370.1, KQ031389.1, KZ559109.1, ML143344.1, ML143350.1, ML143379.1, GL000216.2, KI270879.1, ML143366.1, GL000008.2, GL383520.2, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270830.1, KI270848.1, KI270850.1, KI270862.1, KI270866.1, KQ031384.1, KV575244.1
    ##   - in 'y': GL383578.2, KI270707.1, KI270731.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270733.1, ML143365.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270754.1, KN538361.1, ML143380.1, GL383526.1, KI270895.1, KI270903.1, GL000216.2, KI270879.1, ML143366.1, GL000008.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270848.1, KQ031384.1, KV575244.1, KI270707.1, KI270731.1
    ##   - in 'y': GL877875.1, KI270721.1, KI270761.1, KI270806.1, KI270844.1, KI270845.1, KI270899.1, KI270937.1, KN196479.1, KQ458383.1, KZ208912.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270733.1, ML143365.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270754.1, KN538361.1, ML143380.1, GL383526.1, KI270895.1, KI270903.1, GL000216.2, KI270879.1, ML143366.1, GL000008.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270848.1, KQ031384.1, KV575244.1, KI270707.1, KI270731.1
    ##   - in 'y': GL877875.1, KI270721.1, KI270761.1, KI270806.1, KI270844.1, KI270845.1, KI270899.1, KI270937.1, KN196479.1, KQ458383.1, KZ208912.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270733.1, KI270745.1, KI270818.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270782.1, KN538361.1, KZ559112.1, ML143380.1, GL383526.1, KI270784.1, KI270831.1, KI270880.1, KI270895.1, KI270903.1, KZ559109.1, ML143350.1, GL000216.2, ML143366.1, GL000008.2, GL000224.1, GL000254.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270848.1, KI270862.1, KI270866.1, KQ031384.1, GL383578.2, KI270731.1, KI270721.1, KI270761.1, KI270806.1, KI270844.1, KI270845.1, KI270899.1, KI270937.1, KN196479.1, KQ458383.1, ML143358.1
    ##   - in 'y': KI270816.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000220.1, GL000225.1, GL000251.2, KI270310.1, KI270330.1, KI270438.1, KI270709.1, KI270733.1, KI270745.1, KI270818.1, ML143365.1, ML143371.1, ML143372.1, GL000009.2, GL000214.1, GL000218.1, KI270711.1, KI270782.1, KN538361.1, KZ559112.1, ML143380.1, GL383526.1, KI270784.1, KI270831.1, KI270880.1, KI270895.1, KI270903.1, KZ559109.1, ML143350.1, GL000216.2, ML143366.1, GL000008.2, GL000224.1, GL000254.2, GL383520.2, KI270538.1, KI270729.1, KI270730.1, KI270736.1, KI270751.1, KI270757.1, KI270848.1, KI270862.1, KI270866.1, KQ031384.1, GL383578.2, KI270731.1, KI270721.1, KI270761.1, KI270806.1, KI270844.1, KI270845.1, KI270899.1, KI270937.1, KN196479.1, KQ458383.1, ML143358.1
    ##   - in 'y': KI270816.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270337.1, ML143344.1, ML143352.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000220.1, GL000224.1, GL000256.2, KI270754.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270337.1, ML143344.1, ML143352.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000220.1, GL000224.1, GL000256.2, KI270754.1, KQ031389.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000219.1, GL000251.2, GL000255.2, GL383563.3, GL949752.1, KI270442.1, KI270538.1, KI270714.1, KI270721.1, KI270744.1, KI270816.1, KV766198.1, KZ208915.1, ML143345.1, ML143355.1, ML143366.1
    ##   - in 'y': GL000214.1, GL000224.1, KI270310.1, KI270411.1, KI270438.1, KI270733.1, KI270754.1, KI270849.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000216.2, GL000219.1, GL000251.2, GL000255.2, GL383563.3, GL949752.1, KI270442.1, KI270538.1, KI270714.1, KI270721.1, KI270744.1, KI270816.1, KV766198.1, KZ208915.1, ML143345.1, ML143355.1, ML143366.1
    ##   - in 'y': GL000214.1, GL000224.1, KI270310.1, KI270411.1, KI270438.1, KI270733.1, KI270754.1, KI270849.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270733.1, KI270754.1, KI270936.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL339449.2, GL383533.1, GL383555.2, GL383556.1, GL383578.2, JH159136.1, KI270418.1, KI270435.1, KI270442.1, KI270511.1, KI270519.1, KI270707.1, KI270709.1, KI270710.1, KI270711.1, KI270712.1, KI270714.1, KI270717.1, KI270719.1, KI270720.1, KI270722.1, KI270731.1, KI270732.1, KI270735.1, KI270738.1, KI270742.1, KI270743.1, KI270744.1, KI270746.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270763.1, KI270770.1, KI270772.1, KI270782.1, KI270783.1, KI270787.1, KI270803.1, KI270821.1, KI270822.1, KI270831.1, KI270832.1, KI270836.1, KI270844.1, KI270850.1, KI270854.1, KI270875.1, KI270876.1, KI270878.1, KI270896.1, KI270905.1, KI270908.1, KI270913.1, KI270925.1, KI270927.1, KN196478.1, KN196484.1, KN196487.1, KN538372.1, KQ090016.1, KQ090022.1, KQ090026.1, KV766198.1, KZ208907.1, KZ208908.1, KZ208913.1, KZ208915.1, KZ208922.1, KZ559103.1, KZ559110.1, KZ559112.1, ML143355.1, ML143364.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270733.1, KI270754.1, KI270936.1
    ##   - in 'y': GL000008.2, GL000009.2, GL000194.1, GL000205.2, GL000214.1, GL000216.2, GL000218.1, GL000220.1, GL000221.1, GL000224.1, GL000225.1, GL000253.2, GL000254.2, GL000255.2, GL000256.2, GL000257.2, GL339449.2, GL383533.1, GL383555.2, GL383556.1, GL383578.2, JH159136.1, KI270418.1, KI270435.1, KI270442.1, KI270511.1, KI270519.1, KI270707.1, KI270709.1, KI270710.1, KI270711.1, KI270712.1, KI270714.1, KI270717.1, KI270719.1, KI270720.1, KI270722.1, KI270731.1, KI270732.1, KI270735.1, KI270738.1, KI270742.1, KI270743.1, KI270744.1, KI270746.1, KI270749.1, KI270750.1, KI270751.1, KI270753.1, KI270763.1, KI270770.1, KI270772.1, KI270782.1, KI270783.1, KI270787.1, KI270803.1, KI270821.1, KI270822.1, KI270831.1, KI270832.1, KI270836.1, KI270844.1, KI270850.1, KI270854.1, KI270875.1, KI270876.1, KI270878.1, KI270896.1, KI270905.1, KI270908.1, KI270913.1, KI270925.1, KI270927.1, KN196478.1, KN196484.1, KN196487.1, KN538372.1, KQ090016.1, KQ090022.1, KQ090026.1, KV766198.1, KZ208907.1, KZ208908.1, KZ208913.1, KZ208915.1, KZ208922.1, KZ559103.1, KZ559110.1, KZ559112.1, ML143355.1, ML143364.1, ML143366.1, ML143371.1, ML143372.1, ML143377.1, ML143378.1, chrM
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, GL000008.2, GL000194.1, GL000205.2, GL000221.1, GL000253.2, GL000254.2, GL000257.2, GL339449.2, GL383533.1, GL383555.2, GL383556.1, GL383578.2, JH159136.1, KI270418.1, KI270435.1, KI270511.1, KI270519.1, KI270710.1, KI270717.1, KI270719.1, KI270720.1, KI270722.1, KI270732.1, KI270735.1, KI270738.1, KI270743.1, KI270746.1, KI270749.1, KI270750.1, KI270753.1, KI270763.1, KI270770.1, KI270772.1, KI270782.1, KI270783.1, KI270787.1, KI270822.1, KI270832.1, KI270836.1, KI270844.1, KI270854.1, KI270875.1, KI270876.1, KI270878.1, KI270896.1, KI270913.1, KI270925.1, KI270927.1, KN196478.1, KN196487.1, KN538372.1, KQ090016.1, KQ090022.1, KV766198.1, KZ208907.1, KZ208908.1, KZ208913.1, KZ208922.1, KZ559110.1, KZ559112.1, ML143366.1, ML143378.1, chrM
    ##   - in 'y': GL383545.1, GL949752.1, KI270330.1, KI270730.1, KI270741.1, KI270802.1, KI270804.1, KI270819.1, KI270830.1, KI270848.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270934.1, KN196480.1, KN538370.1, KQ031389.1, KV880763.1, KZ208912.1, ML143345.1, ML143353.1, ML143358.1, ML143375.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, GL000008.2, GL000194.1, GL000205.2, GL000221.1, GL000253.2, GL000254.2, GL000257.2, GL339449.2, GL383533.1, GL383555.2, GL383556.1, GL383578.2, JH159136.1, KI270418.1, KI270435.1, KI270511.1, KI270519.1, KI270710.1, KI270717.1, KI270719.1, KI270720.1, KI270722.1, KI270732.1, KI270735.1, KI270738.1, KI270743.1, KI270746.1, KI270749.1, KI270750.1, KI270753.1, KI270763.1, KI270770.1, KI270772.1, KI270782.1, KI270783.1, KI270787.1, KI270822.1, KI270832.1, KI270836.1, KI270844.1, KI270854.1, KI270875.1, KI270876.1, KI270878.1, KI270896.1, KI270913.1, KI270925.1, KI270927.1, KN196478.1, KN196487.1, KN538372.1, KQ090016.1, KQ090022.1, KV766198.1, KZ208907.1, KZ208908.1, KZ208913.1, KZ208922.1, KZ559110.1, KZ559112.1, ML143366.1, ML143378.1, chrM
    ##   - in 'y': GL383545.1, GL949752.1, KI270330.1, KI270730.1, KI270741.1, KI270802.1, KI270804.1, KI270819.1, KI270830.1, KI270848.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270934.1, KN196480.1, KN538370.1, KQ031389.1, KV880763.1, KZ208912.1, ML143345.1, ML143353.1, ML143358.1, ML143375.1, chrY
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270759.1, KI270936.1, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL000254.2, GL000257.2, GL383555.2, KI270418.1, KI270511.1, KI270709.1, KI270711.1, KI270717.1, KI270722.1, KI270732.1, KI270735.1, KI270738.1, KI270743.1, KI270746.1, KI270749.1, KI270750.1, KI270763.1, KI270770.1, KI270772.1, KI270836.1, KI270925.1, KN196478.1, KN196484.1, KN196487.1, KN538372.1, KQ090026.1, KZ208907.1, KZ208908.1, KZ559110.1, ML143378.1, GL383545.1, GL949752.1, KI270330.1, KI270730.1, KI270741.1, KI270802.1, KI270804.1, KI270819.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270934.1, KN196480.1, KN538370.1, KV880763.1, ML143358.1, ML143375.1, chrY
    ##   - in 'y': GL000250.2, GL000251.2, GL383574.1, KI270718.1, KI270745.1, KI270765.1, KI270779.1, KI270827.1, KI270849.1, KI270851.1, KI270856.1, KI270899.1, KI270904.1, KN538360.1, KQ090023.1, KQ458382.1, KQ983255.1, KV880764.1, ML143344.1, ML143352.1, ML143362.1, ML143365.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270759.1, KI270936.1, GL000216.2, GL000218.1, GL000220.1, GL000225.1, GL000254.2, GL000257.2, GL383555.2, KI270418.1, KI270511.1, KI270709.1, KI270711.1, KI270717.1, KI270722.1, KI270732.1, KI270735.1, KI270738.1, KI270743.1, KI270746.1, KI270749.1, KI270750.1, KI270763.1, KI270770.1, KI270772.1, KI270836.1, KI270925.1, KN196478.1, KN196484.1, KN196487.1, KN538372.1, KQ090026.1, KZ208907.1, KZ208908.1, KZ559110.1, ML143378.1, GL383545.1, GL949752.1, KI270330.1, KI270730.1, KI270741.1, KI270802.1, KI270804.1, KI270819.1, KI270880.1, KI270892.1, KI270897.1, KI270903.1, KI270934.1, KN196480.1, KN538370.1, KV880763.1, ML143358.1, ML143375.1, chrY
    ##   - in 'y': GL000250.2, GL000251.2, GL383574.1, KI270718.1, KI270745.1, KI270765.1, KI270779.1, KI270827.1, KI270849.1, KI270851.1, KI270856.1, KI270899.1, KI270904.1, KN538360.1, KQ090023.1, KQ458382.1, KQ983255.1, KV880764.1, ML143344.1, ML143352.1, ML143362.1, ML143365.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000254.2, GL000255.2, GL000256.2, KI270714.1, KI270728.1, KI270765.1, KI270784.1, KI270805.1, KI270849.1, KQ090026.1, KV575244.1, KV766198.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270709.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000254.2, GL000255.2, GL000256.2, KI270714.1, KI270728.1, KI270765.1, KI270784.1, KI270805.1, KI270849.1, KQ090026.1, KV575244.1, KV766198.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270709.1, KN196487.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, KI270310.1, KN196487.1
    ##   - in 'y': KI270745.1, KI270754.1, KI270831.1, KI270853.1, KI270892.1, KI270905.1, KQ458385.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, KI270310.1, KN196487.1
    ##   - in 'y': KI270745.1, KI270754.1, KI270831.1, KI270853.1, KI270892.1, KI270905.1, KQ458385.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000254.2, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270728.1, KI270733.1, KI270744.1, KI270765.1, KI270784.1, KI270805.1, KI270849.1, KQ090026.1, KV575244.1, KV766198.1, KZ208915.1, chr18, chr21, GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270709.1, KN196487.1, KI270745.1, KI270754.1, KI270831.1, KI270853.1, KI270892.1, KI270905.1, KQ458385.1
    ##   - in 'y': KI270442.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000254.2, GL000255.2, GL000256.2, KI270438.1, KI270714.1, KI270728.1, KI270733.1, KI270744.1, KI270765.1, KI270784.1, KI270805.1, KI270849.1, KQ090026.1, KV575244.1, KV766198.1, KZ208915.1, chr18, chr21, GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, KI270310.1, KI270709.1, KN196487.1, KI270745.1, KI270754.1, KI270831.1, KI270853.1, KI270892.1, KI270905.1, KQ458385.1
    ##   - in 'y': KI270442.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL383563.3, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270718.1, KI270728.1, KI270742.1, KI270751.1, KI270759.1, KI270765.1, KI270782.1, KI270849.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000220.1, KI270438.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL383563.3, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270718.1, KI270728.1, KI270742.1, KI270751.1, KI270759.1, KI270765.1, KI270782.1, KI270849.1, ML143352.1, ML143355.1, ML143379.1
    ##   - in 'y': GL000216.2, GL000220.1, KI270438.1, KI270709.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL339449.2, GL383563.3, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270718.1, KI270751.1, KI270759.1, KI270765.1, KI270849.1, ML143379.1, GL000216.2, GL000220.1, KI270438.1
    ##   - in 'y': KI270899.1, KN538370.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000253.2, GL000255.2, GL339449.2, GL383563.3, KI270337.1, KI270466.1, KI270467.1, KI270712.1, KI270713.1, KI270718.1, KI270751.1, KI270759.1, KI270765.1, KI270849.1, ML143379.1, GL000216.2, GL000220.1, KI270438.1
    ##   - in 'y': KI270899.1, KN538370.1, ML143375.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, GL383526.1, GL383578.2, KI270466.1, KI270467.1, KI270840.1, KI270849.1, KI270860.1, KI270861.1, KV766196.1, KZ559109.1
    ##   - in 'y': GL000214.1, GL949752.1, KI270538.1, KI270726.1, KI270830.1, KI270907.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000254.2, GL383526.1, GL383578.2, KI270466.1, KI270467.1, KI270840.1, KI270849.1, KI270860.1, KI270861.1, KV766196.1, KZ559109.1
    ##   - in 'y': GL000214.1, GL949752.1, KI270538.1, KI270726.1, KI270830.1, KI270907.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270712.1, KI270714.1, KI270718.1, KI270719.1, KI270726.1, KI270728.1, KI270742.1, KI270743.1, KI270754.1, KI270765.1, KI270783.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270855.1, KI270857.1, KI270860.1, KI270908.1, KI270924.1, KN196484.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143344.1, ML143345.1, ML143350.1, ML143352.1, ML143366.1, ML143370.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000220.1, GL000225.1, KI270442.1, KI270589.1, KI270709.1, KI270729.1, KI270736.1, KI270745.1, KI270751.1, KI270805.1, KI270831.1, KI270879.1, KI270880.1, KN196487.1, KN538364.1, KV880768.1, ML143341.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000221.1, GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270712.1, KI270714.1, KI270718.1, KI270719.1, KI270726.1, KI270728.1, KI270742.1, KI270743.1, KI270754.1, KI270765.1, KI270783.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270855.1, KI270857.1, KI270860.1, KI270908.1, KI270924.1, KN196484.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143344.1, ML143345.1, ML143350.1, ML143352.1, ML143366.1, ML143370.1, ML143371.1, ML143377.1
    ##   - in 'y': GL000220.1, GL000225.1, KI270442.1, KI270589.1, KI270709.1, KI270729.1, KI270736.1, KI270745.1, KI270751.1, KI270805.1, KI270831.1, KI270879.1, KI270880.1, KN196487.1, KN538364.1, KV880768.1, ML143341.1, ML143358.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270726.1, KI270728.1, KI270743.1, KI270754.1, KI270765.1, KI270782.1, KI270783.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270855.1, KI270857.1, KI270860.1, KI270894.1, KI270905.1, KI270908.1, KI270924.1, KN196484.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143344.1, ML143345.1, ML143350.1, ML143352.1, ML143365.1, ML143366.1, ML143370.1, ML143371.1, ML143372.1, ML143377.1, ML143379.1, GL000225.1, KI270442.1, KI270589.1, KI270729.1, KI270736.1, KI270751.1, KI270805.1, KI270879.1, KV880768.1, ML143341.1, ML143358.1
    ##   - in 'y': GL000216.2, KI270711.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000008.2, GL000009.2, GL000194.1, GL000214.1, GL000221.1, GL000224.1, GL000251.2, GL000253.2, GL000255.2, GL000256.2, GL383578.2, GL949752.1, KI270706.1, KI270712.1, KI270713.1, KI270714.1, KI270718.1, KI270719.1, KI270726.1, KI270728.1, KI270743.1, KI270754.1, KI270765.1, KI270782.1, KI270783.1, KI270804.1, KI270819.1, KI270821.1, KI270830.1, KI270844.1, KI270849.1, KI270853.1, KI270855.1, KI270857.1, KI270860.1, KI270894.1, KI270905.1, KI270908.1, KI270924.1, KN196484.1, KV766198.1, KV880764.1, KZ208915.1, KZ559112.1, ML143344.1, ML143345.1, ML143350.1, ML143352.1, ML143365.1, ML143366.1, ML143370.1, ML143371.1, ML143372.1, ML143377.1, ML143379.1, GL000225.1, KI270442.1, KI270589.1, KI270729.1, KI270736.1, KI270751.1, KI270805.1, KI270879.1, KV880768.1, ML143341.1, ML143358.1
    ##   - in 'y': GL000216.2, KI270711.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000251.2, GL383578.2, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270712.1, KI270714.1, KI270718.1, KI270733.1, KI270744.1, KI270765.1, KI270783.1, KI270821.1, KI270830.1, KI270849.1, KI270855.1, KI270894.1, KI270924.1, KV766198.1, KZ559112.1, ML143344.1, ML143350.1, ML143365.1, ML143366.1, ML143370.1, ML143372.1, GL000220.1, GL000225.1, KI270442.1, KI270589.1, KI270709.1, KI270729.1, KI270736.1, KI270745.1, KI270751.1, KI270805.1, KI270831.1, KI270879.1, KN196487.1, KN538364.1, KV880768.1, ML143341.1, ML143358.1, GL000216.2, KI270711.1
    ##   - in 'y': KI270734.1, KI270897.1, KI270899.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000224.1, GL000251.2, GL383578.2, KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270712.1, KI270714.1, KI270718.1, KI270733.1, KI270744.1, KI270765.1, KI270783.1, KI270821.1, KI270830.1, KI270849.1, KI270855.1, KI270894.1, KI270924.1, KV766198.1, KZ559112.1, ML143344.1, ML143350.1, ML143365.1, ML143366.1, ML143370.1, ML143372.1, GL000220.1, GL000225.1, KI270442.1, KI270589.1, KI270709.1, KI270729.1, KI270736.1, KI270745.1, KI270751.1, KI270805.1, KI270831.1, KI270879.1, KN196487.1, KN538364.1, KV880768.1, ML143341.1, ML143358.1, GL000216.2, KI270711.1
    ##   - in 'y': KI270734.1, KI270897.1, KI270899.1, KV575244.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270899.1, KI270908.1, ML143369.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270713.1, KI270714.1, KI270733.1, KI270738.1, KI270744.1, KI270751.1, KI270754.1, KI270757.1, KI270821.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL582966.2, KI270337.1, KI270466.1, KI270467.1, KI270899.1, KI270908.1, ML143369.1, ML143380.1
    ##   - in 'y': GL000214.1, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270713.1, KI270714.1, KI270733.1, KI270738.1, KI270744.1, KI270751.1, KI270754.1, KI270757.1, KI270821.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, KI270337.1, KI270466.1, KI270467.1, KI270899.1, KQ031384.1, ML143369.1, ML143380.1, GL000216.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270738.1, KI270751.1, KI270754.1, KI270757.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1
    ##   - in 'y': GL000256.2, KI270442.1, KI270538.1, KI270892.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000194.1, GL000219.1, KI270337.1, KI270466.1, KI270467.1, KI270899.1, KQ031384.1, ML143369.1, ML143380.1, GL000216.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270738.1, KI270751.1, KI270754.1, KI270757.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1
    ##   - in 'y': GL000256.2, KI270442.1, KI270538.1, KI270892.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL582966.2, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270908.1, KQ031384.1, ML143380.1, chr14, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270713.1, KI270714.1, KI270733.1, KI270738.1, KI270744.1, KI270751.1, KI270757.1, KI270821.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1, GL000256.2, KI270442.1, KI270538.1, KI270892.1
    ##   - in 'y': KI270706.1, KI270712.1, KI270728.1, KQ031389.1, KV880764.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000219.1, GL582966.2, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270908.1, KQ031384.1, ML143380.1, chr14, GL000216.2, GL000220.1, GL000224.1, GL000225.1, GL000255.2, GL949752.1, KI270310.1, KI270465.1, KI270707.1, KI270713.1, KI270714.1, KI270733.1, KI270738.1, KI270744.1, KI270751.1, KI270757.1, KI270821.1, KI270868.1, KN196487.1, KZ208915.1, ML143367.1, GL000256.2, KI270442.1, KI270538.1, KI270892.1
    ##   - in 'y': KI270706.1, KI270712.1, KI270728.1, KQ031389.1, KV880764.1, ML143355.1, ML143379.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270728.1, KI270735.1
    ##   - in 'y': GL000009.2, GL000218.1, GL000219.1, GL000220.1, GL000250.2, GL000251.2, GL000253.2, GL383579.2, KI270706.1, KI270717.1, KI270726.1, KI270733.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270772.1, KI270799.1, KI270830.1, KI270831.1, KI270842.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270905.1, KI270908.1, KN196487.1, KQ458383.1, KV766193.1, KZ559105.1, ML143345.1, ML143350.1, ML143366.1, ML143367.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270466.1, KI270467.1, KI270709.1, KI270728.1, KI270735.1
    ##   - in 'y': GL000009.2, GL000218.1, GL000219.1, GL000220.1, GL000250.2, GL000251.2, GL000253.2, GL383579.2, KI270706.1, KI270717.1, KI270726.1, KI270733.1, KI270742.1, KI270743.1, KI270744.1, KI270745.1, KI270772.1, KI270799.1, KI270830.1, KI270831.1, KI270842.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270905.1, KI270908.1, KN196487.1, KQ458383.1, KV766193.1, KZ559105.1, ML143345.1, ML143350.1, ML143366.1, ML143367.1, ML143371.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270337.1, KI270466.1, KI270467.1, KI270728.1, KI270735.1, GL000009.2, GL000218.1, GL000220.1, GL000250.2, KI270706.1, KI270717.1, KI270726.1, KI270742.1, KI270743.1, KI270745.1, KI270799.1, KI270830.1, KI270831.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270908.1, KQ458383.1, KV766193.1, ML143345.1, ML143350.1, ML143371.1
    ##   - in 'y': GL000216.2, GL000224.1, GL000225.1, KI270310.1, KI270442.1, KI270538.1, KI270729.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000256.2, KI270337.1, KI270466.1, KI270467.1, KI270728.1, KI270735.1, GL000009.2, GL000218.1, GL000220.1, GL000250.2, KI270706.1, KI270717.1, KI270726.1, KI270742.1, KI270743.1, KI270745.1, KI270799.1, KI270830.1, KI270831.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270880.1, KI270908.1, KQ458383.1, KV766193.1, ML143345.1, ML143350.1, ML143371.1
    ##   - in 'y': GL000216.2, GL000224.1, GL000225.1, KI270310.1, KI270442.1, KI270538.1, KI270729.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270719.1, KI270728.1, KI270735.1, ML143380.1, GL000009.2, GL000218.1, GL000219.1, GL000220.1, GL000250.2, GL000253.2, GL383579.2, KI270706.1, KI270717.1, KI270726.1, KI270743.1, KI270745.1, KI270772.1, KI270799.1, KI270830.1, KI270831.1, KI270842.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270905.1, KI270908.1, KN196487.1, KQ458383.1, KV766193.1, KZ559105.1, ML143345.1, ML143350.1, ML143366.1, ML143367.1, ML143371.1, GL000216.2, GL000224.1, GL000225.1, KI270310.1, KI270442.1, KI270538.1, KI270729.1
    ##   - in 'y': GL949752.1, KI270714.1, KI270754.1, KI270819.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': KI270337.1, KI270438.1, KI270466.1, KI270467.1, KI270709.1, KI270712.1, KI270719.1, KI270728.1, KI270735.1, ML143380.1, GL000009.2, GL000218.1, GL000219.1, GL000220.1, GL000250.2, GL000253.2, GL383579.2, KI270706.1, KI270717.1, KI270726.1, KI270743.1, KI270745.1, KI270772.1, KI270799.1, KI270830.1, KI270831.1, KI270842.1, KI270851.1, KI270853.1, KI270856.1, KI270878.1, KI270905.1, KI270908.1, KN196487.1, KQ458383.1, KV766193.1, KZ559105.1, ML143345.1, ML143350.1, ML143366.1, ML143367.1, ML143371.1, GL000216.2, GL000224.1, GL000225.1, KI270310.1, KI270442.1, KI270538.1, KI270729.1
    ##   - in 'y': GL949752.1, KI270714.1, KI270754.1, KI270819.1, ML143355.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270438.1, KI270466.1, KI270467.1, KN538372.1, ML143379.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, KI270731.1, KI270733.1, KI270742.1, KI270754.1, KI270853.1, KI270879.1, KI270903.1, KZ208915.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

    ## Warning in .Seqinfo.mergexy(x, y): Each of the 2 combined objects has sequence levels not in the other:
    ##   - in 'x': GL000255.2, KI270438.1, KI270466.1, KI270467.1, KN538372.1, ML143379.1
    ##   - in 'y': GL000009.2, GL000194.1, GL000214.1, GL000216.2, GL000218.1, GL000219.1, GL000224.1, KI270731.1, KI270733.1, KI270742.1, KI270754.1, KI270853.1, KI270879.1, KI270903.1, KZ208915.1, KZ559109.1
    ##   Make sure to always combine/compare objects based on the same reference
    ##   genome (use suppressWarnings() to suppress this warning).

Number of peaks threshold
-------------------------

Since we want to study robust DNA binding proteins in K562 cells, we will implement a cutoff for the minimum number of peaks a DBP must have in order to be considered in the following analyses. There are some DBPs have have no replicate concordant peaks.

``` r
num_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))
g <- ggplot(num_peaks_df, aes(x = num_peaks))
g + geom_histogram(bins = 70) +
  xlab("Number of consensus peaks") +
  ylab("Count") +
  ggtitle("Distribution of number of consensus peaks")
```

![](consensus_peak_set_files/figure-markdown_github/dist-num-peaks-1.png)

``` r
ggsave("figures/consensus_peaks_histogram.pdf")
```

    ## Saving 7 x 5 in image

#### No replicate concordance

The proteins with zero consensus peaks are: `cat(paste(num_peaks_df[which(num_peaks_df$num_peaks == 0), "dbp"], collapse = " "))`

#### Peaks cutoff

``` r
# We're going toconsensus_peaks apply a cutoff at 250 peaks
num_peaks_threshold <- 250
consensus_peaks <- consensus_peaks[num_peaks_df$num_peaks > num_peaks_threshold]
```

Since this captures the majority of DPBs and still provides a reasonable number of peaks to work with, we chose a cutoff of `num_peaks_threshold` peaks. This results in losing the following proteins: `cat(paste(num_peaks_df[which(num_peaks_df$num_peaks <= 250), "dbp"], collapse = " "))`

``` r
# Export the peak lists.
for(i in 1:length(consensus_peaks)) {
  rtracklayer::export(consensus_peaks[[i]], paste0("results/", names(consensus_peaks)[i], "_consensus_peaks.bed"))
}
```

Summary of consensus peaks
--------------------------

Here we'll look at a few characteristics of the remaining peaksets.

#### Total peak length

``` r
# Subset to remaining peaks
num_peaks_df <- num_peaks_df %>% filter(dbp %in% names(consensus_peaks))

# Calculate the total peak width (bp bound by all peaks)
num_peaks_df$total_peak_length <- sapply(consensus_peaks, function(peaks) sum(width(peaks)))

g <- ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length, label = dbp))
g + geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "black", lty = 2) +
  geom_text() +
  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](consensus_peak_set_files/figure-markdown_github/total-peak-length-1.png)

``` r
ggsave("figures/peak_count_vs_peak_length.pdf")
```

    ## Saving 7 x 5 in image

#### Peak width distributions

``` r
# Let's grab the peak widths for all peak sets
peak_widths_df <- lapply(consensus_peaks, 
                      function(peaks) paste(width(peaks), collapse = ";")) %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = "dbp",
               values_to = "peak_width") %>%
  separate_rows("peak_width", sep = ";", convert = T)

peak_widths_summary <- peak_widths_df %>% 
  group_by(dbp) %>%
  summarize("mean_width" = mean(peak_width),
            "median_width" = median(peak_width),
            "sd_width" = sd(peak_width),
            "cv_width" = (sd_width/mean_width)*100)

g <- ggplot(peak_widths_summary, aes(x = mean_width))
g + geom_histogram(bins = 60) + 
  xlab("Mean peak width") +
  ylab("Count") +
  ggtitle("Peak width distribution")
```

![](consensus_peak_set_files/figure-markdown_github/peak-width-distribution-1.png)

``` r
ggsave("figures/peak_width_distribution.pdf")
```

    ## Saving 7 x 5 in image

``` r
g <- ggplot(peak_widths_summary, aes(x = mean_width, y = cv_width, label = dbp))
g + geom_point() + 
  geom_text() + 
  xlab("Mean peak width") +
  ylab("CV peak width") +
  ggtitle("Peak width: mean vs. coefficient of variation")
```

![](consensus_peak_set_files/figure-markdown_github/peak-width-distribution-2.png)

``` r
ggsave("figures/peak_width_vs_cv.pdf")
```

    ## Saving 7 x 5 in image

It seems that most DBPs have a relatively small peak width ~ 1000bp +/- 500bp and relatively small CV ~100bp. However, a few including POLII subunits have wider peaks widths. RFX1 is unique in that it has a typical mean peak width, but a very large CV. Indicating that while most peaks are in the typical range, it may have some very large peaks.

``` r
# Let's look at PolII in particular
g <- ggplot(peak_widths_df %>% filter(dbp %in% c("POLR2A", "POLR2B", "SUPT5H")), aes(x = log10(peak_width)))
g + geom_histogram(bins = 100) + 
  facet_grid(dbp~., scales = "free_y") + 
  xlab("log10(Peak width)") +
  ylab("Count") +
  ggtitle("Peak width distribution: POLII")
```

![](consensus_peak_set_files/figure-markdown_github/polII-peak-widths-1.png)

``` r
ggsave("figures/peak_width_distribution_polII.pdf")
```

    ## Saving 7 x 5 in image

So interestingly it seems that there are two peaks in POLR2A's binding, but not for the other subunits. This may indicate two binding modes for POLR2A.

``` r
# Also RFX1 which had a very high CV.
g <- ggplot(peak_widths_df %>% filter(dbp %in% c("RFX1")), aes(x = log10(peak_width)))
g + geom_histogram(bins = 100) + 
  facet_grid(dbp~., scales = "free_y") + 
  xlab("log10(Peak width)") +
  ylab("Count") +
  ggtitle("Peak width distribution: RFX1")
```

![](consensus_peak_set_files/figure-markdown_github/rfx1-peak-widths-1.png)

``` r
ggsave("figures/peak_width_distribution_RFX1.pdf")
```

    ## Saving 7 x 5 in image

``` r
max(peak_widths_df[which(peak_widths_df$dbp == "RFX1" & peak_widths_df$peak_width > 3000),"peak_width"])
```

    ## [1] 286683

And it does seem that RFX1 has just `nrow(peak_widths_df[which(peak_widths_df$dbp == "RFX1" & peak_widths_df$peak_width > 3000),])` peaks which are above 3000 bps which is skewing the distribution. With one peak that has a width of `max(peak_widths_df[which(peak_widths_df$dbp == "RFX1" & peak_widths_df$peak_width > 3000),"peak_width"])`.

Cutdapt Light Trap OTUs300
================
Hendra Sihaloho

- [Import otu table](#import-otu-table)
- [Using Shell Terminal](#using-shell-terminal)
- [Meta](#meta)
- [Refseq](#refseq)
- [Phyloseq](#phyloseq)

``` r
suppressPackageStartupMessages({
  library(tidyverse) ; library(reshape2) ; library(Biostrings) ; library(phyloseq)
})
```

    ## Warning: package 'ggplot2' was built under R version 4.3.3

# Import otu table

``` r
otu = read.table("./output/otus_finalize/otutab.300.txt", row.names = 1)
length(which(colSums(otu) < 10000))/ncol(otu)
```

    ## [1] 0

``` r
hist(colSums(otu))
```

![](otus_gfm_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
head(otu)
```

    ##          V2   V3   V4  V5  V6  V7   V8   V9  V10 V11 V12 V13 V14  V15  V16  V17
    ## Otu70   222    0    0   0   0   0    0  240    0   0   0   0   0    0  368    1
    ## Otu3   2336 1810 2066 313 701 313 1225 3197 2962 514 313 827 192 1080 2060 2717
    ## Otu19   115  254  277  14 242 108  121  134  406  32   4 155 102  172  101  489
    ## Otu245   40    0    1   0   0   0    0   66    0   0   0   0   0    0   96    0
    ## Otu77    47   44  211   0  65   0   63  257  165   0   0  49   0  102   32  229
    ## Otu103   83    0    0   0   3   0   30   69    0   0   0   3   0   56   71    0
    ##        V18 V19 V20 V21  V22  V23  V24 V25 V26  V27 V28  V29  V30  V31 V32  V33
    ## Otu70    0   0   0   0    0  293    0   0   0    0   0    3  173    0   0    0
    ## Otu3   386 322 955 145 1512 3422 1681 464 837 1331 180 1335 2247 2745 526 1576
    ## Otu19    8  15 153  87  222  263  338  44 179  171  58  257  146  146   3  164
    ## Otu245   0   0   0   0    0   78    0   0   0    0   0    0   44    0   0    0
    ## Otu77    0   0 172   0   52   35  116   0  75   72   0   38   29  256   0  114
    ## Otu103   0   0   3   0   65  177    0   0   5    4   0  101   53    1   0    6
    ##        V34 V35  V36  V37  V38 V39 V40 V41 V42  V43  V44  V45 V46 V47 V48 V49
    ## Otu70    0   0    0  324    0   0   0   0   0    0  189    0   0   0   0   0
    ## Otu3   572 156 1306 2862 2397 927 892 136 146 1583 2752 1915 363 950 162 172
    ## Otu19  149  87  236  197  245  10 207 126 120  278  209  328  14 361 182 151
    ## Otu245   0   0    0   92    0   0   0   0   0    0   45    0   0   0   0   0
    ## Otu77  116   0   52   57   46   0  67   0   0  104   76   74   0 480   0   0
    ## Otu103   2   0   69   67    0   0   1   0   0   57   90    0   0   4   0   0
    ##         V50  V51  V52  V53 V54 V55 V56  V57  V58  V59  V60  V61 V62 V63 V64
    ## Otu70     0  435  403    0   0   0   0    0    0  203  213    0   0   0   0
    ## Otu3   2093 3105 2514 1896 301 846 186 1107 1123 2674 2569 1910 478 731 298
    ## Otu19   191  289   87  314  19 257 228  249  247  108  187  408   1 160 110
    ## Otu245    0  337   39    0   0   0   0    0    0   95  115    0   0   0   0
    ## Otu77    42   37   30  153   0  85   0  234   50   41   60  105   0 157   0
    ## Otu103   54  239   75    0   0   5   0  123   22  109   82    0   0   3   0
    ##         V65  V66  V67  V68  V69 V70 V71 V72  V73
    ## Otu70     0    0  300    0    0   0   0   0    0
    ## Otu3   1875 1397 2395 1809 2221 364 788 169 1242
    ## Otu19   225  302  160  259  283  15 215  97  216
    ## Otu245    0    0   45    0    0   0   0   0    0
    ## Otu77   113   46   34   98   88   0 118   0   61
    ## Otu103  160   65   68    0    0   0   5   0   52

``` r
dim(otu)
```

    ## [1] 3247   72

# Using Shell Terminal

``` bash
grep "^#" ./output/otus_finalize/otutab.300.txt > temp_colnames
sed -E 's/#OTU ID\t//' temp_colnames > temp_colnames2
### Using Shell Terminal ===
```

``` r
temp_colnames = read.table("temp_colnames2")
dim(temp_colnames)
```

    ## [1]  1 72

``` r
head(temp_colnames)
```

    ##               V1             V2             V3             V4             V5
    ## 1 COI-TK05D_S460 COI-TK06B_S470 COI-TK06L_S480 COI-TK07J_S490 COI-TK08H_S500
    ##               V6             V7             V8             V9            V10
    ## 1 COI-TK09F_S510 COI-TK10D_S520 COI-TK05E_S461 COI-TK06C_S471 COI-TK07A_S481
    ##              V11            V12            V13            V14            V15
    ## 1 COI-TK07K_S491 COI-TK08I_S501 COI-TK09G_S511 COI-TK10E_S521 COI-TK05F_S462
    ##              V16            V17            V18            V19            V20
    ## 1 COI-TK06D_S472 COI-TK07B_S482 COI-TK07L_S492 COI-TK08J_S502 COI-TK09H_S512
    ##              V21            V22            V23            V24            V25
    ## 1 COI-TK10F_S522 COI-TK05G_S463 COI-TK06E_S473 COI-TK07C_S483 COI-TK08A_S493
    ##              V26            V27            V28            V29            V30
    ## 1 COI-TK08K_S503 COI-TK09I_S513 COI-TK10G_S523 COI-TK05H_S464 COI-TK06F_S474
    ##              V31            V32            V33            V34            V35
    ## 1 COI-TK07D_S484 COI-TK08B_S494 COI-TK08L_S504 COI-TK09J_S514 COI-TK10H_S524
    ##              V36            V37            V38            V39            V40
    ## 1 COI-TK05I_S465 COI-TK06G_S475 COI-TK07E_S485 COI-TK08C_S495 COI-TK09A_S505
    ##              V41            V42            V43            V44            V45
    ## 1 COI-TK09K_S515 COI-TK10I_S525 COI-TK05J_S466 COI-TK06H_S476 COI-TK07F_S486
    ##              V46            V47            V48            V49            V50
    ## 1 COI-TK08D_S496 COI-TK09B_S506 COI-TK09L_S516 COI-TK10J_S526 COI-TK05A_S457
    ##              V51            V52            V53            V54            V55
    ## 1 COI-TK05K_S467 COI-TK06I_S477 COI-TK07G_S487 COI-TK08E_S497 COI-TK09C_S507
    ##              V56            V57            V58            V59            V60
    ## 1 COI-TK10A_S517 COI-TK10K_S527 COI-TK05B_S458 COI-TK05L_S468 COI-TK06J_S478
    ##              V61            V62            V63            V64            V65
    ## 1 COI-TK07H_S488 COI-TK08F_S498 COI-TK09D_S508 COI-TK10B_S518 COI-TK10L_S528
    ##              V66            V67            V68            V69            V70
    ## 1 COI-TK05C_S459 COI-TK06A_S469 COI-TK06K_S479 COI-TK07I_S489 COI-TK08G_S499
    ##              V71            V72
    ## 1 COI-TK09E_S509 COI-TK10C_S519

``` r
temp_colnames = gsub("COI-", "", temp_colnames)
temp_colnames = gsub("_S.*", "", temp_colnames)
temp_colnames = data.frame(t(temp_colnames))

length(colnames(otu))
```

    ## [1] 72

``` r
head(temp_colnames)
```

    ##      X1    X2    X3    X4    X5    X6    X7    X8    X9   X10   X11   X12   X13
    ## 1 TK05D TK06B TK06L TK07J TK08H TK09F TK10D TK05E TK06C TK07A TK07K TK08I TK09G
    ##     X14   X15   X16   X17   X18   X19   X20   X21   X22   X23   X24   X25   X26
    ## 1 TK10E TK05F TK06D TK07B TK07L TK08J TK09H TK10F TK05G TK06E TK07C TK08A TK08K
    ##     X27   X28   X29   X30   X31   X32   X33   X34   X35   X36   X37   X38   X39
    ## 1 TK09I TK10G TK05H TK06F TK07D TK08B TK08L TK09J TK10H TK05I TK06G TK07E TK08C
    ##     X40   X41   X42   X43   X44   X45   X46   X47   X48   X49   X50   X51   X52
    ## 1 TK09A TK09K TK10I TK05J TK06H TK07F TK08D TK09B TK09L TK10J TK05A TK05K TK06I
    ##     X53   X54   X55   X56   X57   X58   X59   X60   X61   X62   X63   X64   X65
    ## 1 TK07G TK08E TK09C TK10A TK10K TK05B TK05L TK06J TK07H TK08F TK09D TK10B TK10L
    ##     X66   X67   X68   X69   X70   X71   X72
    ## 1 TK05C TK06A TK06K TK07I TK08G TK09E TK10C

``` r
colnames(otu)= temp_colnames[1,]

head(otu)
```

    ##        TK05D TK06B TK06L TK07J TK08H TK09F TK10D TK05E TK06C TK07A TK07K TK08I
    ## Otu70    222     0     0     0     0     0     0   240     0     0     0     0
    ## Otu3    2336  1810  2066   313   701   313  1225  3197  2962   514   313   827
    ## Otu19    115   254   277    14   242   108   121   134   406    32     4   155
    ## Otu245    40     0     1     0     0     0     0    66     0     0     0     0
    ## Otu77     47    44   211     0    65     0    63   257   165     0     0    49
    ## Otu103    83     0     0     0     3     0    30    69     0     0     0     3
    ##        TK09G TK10E TK05F TK06D TK07B TK07L TK08J TK09H TK10F TK05G TK06E TK07C
    ## Otu70      0     0   368     1     0     0     0     0     0   293     0     0
    ## Otu3     192  1080  2060  2717   386   322   955   145  1512  3422  1681   464
    ## Otu19    102   172   101   489     8    15   153    87   222   263   338    44
    ## Otu245     0     0    96     0     0     0     0     0     0    78     0     0
    ## Otu77      0   102    32   229     0     0   172     0    52    35   116     0
    ## Otu103     0    56    71     0     0     0     3     0    65   177     0     0
    ##        TK08A TK08K TK09I TK10G TK05H TK06F TK07D TK08B TK08L TK09J TK10H TK05I
    ## Otu70      0     0     0     3   173     0     0     0     0     0     0   324
    ## Otu3     837  1331   180  1335  2247  2745   526  1576   572   156  1306  2862
    ## Otu19    179   171    58   257   146   146     3   164   149    87   236   197
    ## Otu245     0     0     0     0    44     0     0     0     0     0     0    92
    ## Otu77     75    72     0    38    29   256     0   114   116     0    52    57
    ## Otu103     5     4     0   101    53     1     0     6     2     0    69    67
    ##        TK06G TK07E TK08C TK09A TK09K TK10I TK05J TK06H TK07F TK08D TK09B TK09L
    ## Otu70      0     0     0     0     0     0   189     0     0     0     0     0
    ## Otu3    2397   927   892   136   146  1583  2752  1915   363   950   162   172
    ## Otu19    245    10   207   126   120   278   209   328    14   361   182   151
    ## Otu245     0     0     0     0     0     0    45     0     0     0     0     0
    ## Otu77     46     0    67     0     0   104    76    74     0   480     0     0
    ## Otu103     0     0     1     0     0    57    90     0     0     4     0     0
    ##        TK10J TK05A TK05K TK06I TK07G TK08E TK09C TK10A TK10K TK05B TK05L TK06J
    ## Otu70      0   435   403     0     0     0     0     0     0   203   213     0
    ## Otu3    2093  3105  2514  1896   301   846   186  1107  1123  2674  2569  1910
    ## Otu19    191   289    87   314    19   257   228   249   247   108   187   408
    ## Otu245     0   337    39     0     0     0     0     0     0    95   115     0
    ## Otu77     42    37    30   153     0    85     0   234    50    41    60   105
    ## Otu103    54   239    75     0     0     5     0   123    22   109    82     0
    ##        TK07H TK08F TK09D TK10B TK10L TK05C TK06A TK06K TK07I TK08G TK09E TK10C
    ## Otu70      0     0     0     0     0   300     0     0     0     0     0     0
    ## Otu3     478   731   298  1875  1397  2395  1809  2221   364   788   169  1242
    ## Otu19      1   160   110   225   302   160   259   283    15   215    97   216
    ## Otu245     0     0     0     0     0    45     0     0     0     0     0     0
    ## Otu77      0   157     0   113    46    34    98    88     0   118     0    61
    ## Otu103     0     3     0   160    65    68     0     0     0     5     0    52

# Meta

``` r
temp = c("TK05", "TK06", "TK07", "TK08", "TK09", "TK10")
meta = data.frame("TK_number"=colnames(otu))
meta$lt_code = meta$TK_number
meta$lt_code = str_sub(meta$lt_code, 1, 4)
unique(meta$lt_code)
```

    ## [1] "TK05" "TK06" "TK07" "TK08" "TK09" "TK10"

``` r
rownames(meta) = meta$TK_number
head(meta)
```

    ##       TK_number lt_code
    ## TK05D     TK05D    TK05
    ## TK06B     TK06B    TK06
    ## TK06L     TK06L    TK06
    ## TK07J     TK07J    TK07
    ## TK08H     TK08H    TK08
    ## TK09F     TK09F    TK09

# Refseq

``` r
ref_seq = Biostrings::readDNAStringSet("./output/otus_finalize/otus97.final.300.fa")
ref_seq
```

    ## DNAStringSet object of length 3247:
    ##        width seq                                            names               
    ##    [1]   313 ACTTTCATCTAATATTGCCCAT...ATTCTTTATCAACATTTATTT Otu1
    ##    [2]   313 CTTATCTTCCAATATCGCACAT...ATTCTTTATCAACATTTATTT Otu2
    ##    [3]   313 CTTATCTTCTAATATTGCACAC...ATTCTTTATCAACACTTATTT Otu3
    ##    [4]   313 ATTATCAAATAATTTATATCAT...ATTTTATATCAACATTTATTT Otu4
    ##    [5]   313 TTTAGCTTCTAATTCTTTCCAT...ATTTTATATCAACATTTATTT Otu5
    ##    ...   ... ...
    ## [3243]   313 ACTCTCTTCAGGGATTGCTCAT...ATTCTTTATCAACATTTATTC Otu3306
    ## [3244]   313 ATTATCTTCGGGAATTGCTCAT...ATCCTATATCAACATTTATTT Otu3307
    ## [3245]   313 TTTAGCATCTAATACCTTTCAC...ATTCTATACCAACACTTATTT Otu3308
    ## [3246]   313 ACTTTCTTCTAATATTGCTCAT...ATTCTATATCAACATTTATTT Otu3309
    ## [3247]   313 CCTATCAAATAATATATTCCAC...ATTCTTTATCAACACCTATTT Otu3310

``` r
range(ref_seq@ranges@width)
```

    ## [1] 300 313

``` r
ref_seq2 = data.frame(ref_seq)
ref_seq2$length = ref_seq@ranges@width
```

# Phyloseq

``` r
# We need to convert taxa_blast and otu to phyloseq objects to phyloseq object (matrix) before merging
otu2 = otu_table(as.matrix(otu), taxa_are_rows = T)

dim(meta)
```

    ## [1] 72  2

``` r
dim(otu2)
```

    ## [1] 3247   72

``` r
length(ref_seq@ranges@group)
```

    ## [1] 3247

``` r
otu_phy = merge_phyloseq(otu2, ref_seq, sample_data(meta))
otu_phy
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3247 taxa and 72 samples ]
    ## sample_data() Sample Data:       [ 72 samples by 2 sample variables ]
    ## refseq()      DNAStringSet:      [ 3247 reference sequences ]

``` r
# check if otu is correct
length(intersect(rownames(otu_table(otu_phy)), rownames(otu)))
```

    ## [1] 3247

``` r
sessioninfo::session_info(pkgs="attached")
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.3.2 (2023-10-31 ucrt)
    ##  os       Windows 11 x64 (build 22631)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/Chicago
    ##  date     2024-08-22
    ##  pandoc   3.1.11 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version date (UTC) lib source
    ##  BiocGenerics * 0.46.0  2023-04-25 [1] Bioconductor
    ##  Biostrings   * 2.68.1  2023-05-16 [1] Bioconductor
    ##  dplyr        * 1.1.3   2023-09-03 [1] CRAN (R 4.3.2)
    ##  forcats      * 1.0.0   2023-01-29 [1] CRAN (R 4.3.0)
    ##  GenomeInfoDb * 1.36.4  2023-10-02 [1] Bioconductor
    ##  ggplot2      * 3.4.4   2023-10-12 [1] CRAN (R 4.3.3)
    ##  IRanges      * 2.34.1  2023-06-22 [1] Bioconductor
    ##  lubridate    * 1.9.3   2023-09-27 [1] CRAN (R 4.3.1)
    ##  phyloseq     * 1.44.0  2023-04-25 [1] Bioconductor
    ##  purrr        * 1.0.2   2023-08-10 [1] CRAN (R 4.3.1)
    ##  readr        * 2.1.4   2023-02-10 [1] CRAN (R 4.3.2)
    ##  reshape2     * 1.4.4   2020-04-09 [1] CRAN (R 4.3.0)
    ##  S4Vectors    * 0.38.2  2023-09-22 [1] Bioconductor
    ##  stringr      * 1.5.0   2022-12-02 [1] CRAN (R 4.3.2)
    ##  tibble       * 3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
    ##  tidyr        * 1.3.0   2023-01-24 [1] CRAN (R 4.3.2)
    ##  tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.3.2)
    ##  XVector      * 0.40.0  2023-04-25 [1] Bioconductor
    ## 
    ##  [1] C:/Users/Hendra Sihaloho/AppData/Local/R/win-library/4.3
    ##  [2] C:/Program Files/R/R-4.3.2/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────

# Duplex basecalling

Observations, Rates, Sample overview

## Pipeline

After much work and bugfixing, [duplex pipeline](https://github.com/imgag/NCCT_scripts/duplex) now runs relatively stable on SRV025. It is now fully included into the T2T-Pipeline

## Duplex rates

We want to find out the duplex rates, here is a list of "normal" WGS Samples, to analyse:

```txt
23038LRa533_01
23038LRa534_01
23038LRa535_01
23038LRa536_01
23038LRa538_01
23038LRa539_01
23038LRa540_01
23038LRa543_01
23038LRa544_01
23038LRa545_01
21073LRa257_01
21073LRa270_01
21073LRa243_01
21073LRa276_01
21073LRa228_01
21073LRa273_01
21073LRa245L2_01
21073LRa256L2_01
21073LRa258L2_01
21073LRa259L2_01At the same time we
```

After checking the assembly results on data sequenced by ourselvs we see that Duplex does not really improve assembly quality over Herro corrected reads. We decided to not pursue Duplex calling for the real T2T genomes.

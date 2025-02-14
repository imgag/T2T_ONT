---
marp: true
paginate: true
theme: gaia
title: 2024-02-15 IMAGAG longread meeting
--- 

# IMGAG longread meeting

17.12.2024
2nd status update T2T

---

![bg h:700](../img/minknow_metrics_B3.png)

---

![w:550](../img/pores_PBA75688.png) Batch2, 96Gb
![w:550](../img/pores_PBC48743.png) Batch3, T202, 80Gb 
![w:550](../img/pores_PBC38194.png) Batch3, T201, 12Gb
Pores getting blocked very quickly (grey)

---

*Possible explanations*:
- Contamination (Salts, ...)
- Heavy Overloading 

*Unlikely reasons*:
- Underloading
- DNA problems

The three low performing samples were isolated on the same day. 

---
Herro vs published

![bg h:700](../img/herro_vs_duplex_published.png)

--- 

### Our first genome

![bg left:60% h:600](../../assembly/qc/unphased_verkko/TUE_02/bandage_graph.no_colors.png)

- 50x UL
- ~20x Duplex + 20x Herro
- Error rate 0.05%
- NG50: 27Mb
- Rcov: 96.18%
- Rdup: 91.58%
- #breaks: 24924
- T2T haplotypes: 0

---

## Other tests: Use verkko polishing step?

- No impact

--- 

## Next steps:

- Annotate assembly graph
- Add methylation
- Further try out assembly with herro reads (Filter length)
- Try novel assembly approach (A. DiGenova)

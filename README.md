# SRMR

Spatially and Robustly Hybrid Mixture Regression Model for Interence of Spatial Dependence

## Abstract

In this paper, we propose a Spatial Robust Mixture Regression model to investigate the relationship between a response variable and a set of explanatory variables over the spatial domain, assuming that the relationships may exhibit complex spatially dynamic patterns that cannot be captured by constant regression coefficients. Our method integrates the robust finite mixture Gaussian regression model with spatial constraints, to simultaneously handle the spatial nonstationarity, local homogeneity, and outlier contaminations. Compared with existing spatial regression models, our proposed model assumes the existence a few distinct regression models that are estimated based on observations that exhibit similar response-predictor relationships. As such, the proposed model not only accounts for nonstationarity in the spatial trend, but also clusters observations into a few distinct and homogenous groups. This provides an advantage on interpretation with a few stationary sub-processes identified that capture the predominant relationships between response and predictor variables. Moreover, the proposed method incorporates robust procedures to handle contaminations from both regression outliers and spatial outliers. By doing so, we robustly segment the spatial domain into distinct local regions with similar regression coefficients, and sporadic locations that are purely outliers. Rigorous statistical hypothesis testing procedure has been designed to test the significance of such segmentation. Experimental results on many synthetic and real-world datasets demonstrate the robustness, accuracy, and effectiveness of our proposed method, compared with other robust finite mixture regression, spatial regression and spatial segmentation methods.

Our paper is available at arXiv [paper link](https://arxiv.org/abs/2109.00539).

## Figures

Figure 1. Simulation experiments

<img src="https://github.com/changwn/SRMR/blob/main/fig1.png" width="800">

Figure 2. Real world data applications

<img src="https://github.com/changwn/SRMR/blob/main/fig2.png" width="900">

## Citations

If you find the code helpful in your resarch or work, please cite the following papers.
```BibTex
@article{chang2021spatially,
  title={Spatially and Robustly Hybrid Mixture Regression Model for Inference of Spatial Dependence},
  author={Chang, Wennan and Dang, Pengtao and Wan, Changlin and Fang, Yue and Zhao, Tong and Zang, Yong and Li, Bo and Zhang, Chi and Cao, Sha},
  journal={arXiv preprint arXiv:2109.00539},
  year={2021}
}
```

## Questions & Problems

- [Wennan Chang](https://changwn.github.io/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Chi Zhang](https://medicine.iu.edu/departments/genetics/faculty/27057/zhang-chi/)
(czhang87@iu.edu)

Assistant Professor

Department of Medical & Molecular Genetics, Indiana University School of Medicine

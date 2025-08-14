# **Troubleshooting**

Here are some common problems and their solutions. If you cannot find your problem here, please file an issue.

------------------------------------------------------------------------

> [!ERROR]\
> ## Error in normalize.quantiles(dataset0) : ERROR; return code from pthread_create() is 22

To solve this problem, you need to install `preprocessCore` without threading support. Try:

```{shell}
git clone https://github.com/bmbolstad/preprocessCore.git
cd preprocessCore
R CMD INSTALL --configure-args="--disable-threading"  .
```

or

```{r}
BiocManager::install(
  "preprocessCore",
  configure.args = "--disable-threading",
  force = TRUE
)
```

See [bioconductor_docker/issues/22](https://github.com/Bioconductor/bioconductor_docker/issues/22), [Scissor/issues/15](https://github.com/sunduanchen/Scissor/issues/15) for more details.

------------------------------------------------------------------------

> ##
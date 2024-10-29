# LDSC (LD Score Regression) `v2.0.0`

LDSC (LD Score Regression) is a command-line tool for estimating heritability and genetic correlation from GWAS summary statistics. It also computes LD Scores. This tool is essential for researchers in genetics and genomics aiming to understand the genetic architecture of complex traits.

## Table of Contents

- [Background](#background)
- [Scientific Foundation](#scientific-foundation)
- [Installation](#installation)
- [Running LDSC](#running-ldsc)
- [Testing](#testing)
- [Contributing](#contributing)
- [Citing LDSC](#citing-ldsc)
- [License](#license)
- [Authors](#authors)

## Background

Genome-wide association studies (GWAS) have identified thousands of genetic variants associated with complex traits. However, interpreting these associations requires robust statistical tools. LDSC provides a framework to estimate the heritability of traits and the genetic correlation between them using summary statistics from GWAS, leveraging linkage disequilibrium (LD) patterns.

## Scientific Foundation

LDSC implements LD Score regression, a method that distinguishes confounding biases from true polygenic signals in GWAS data. By regressing GWAS test statistics on LD Scores, LDSC estimates the proportion of variance in a trait explained by genetic factors (heritability) and assesses the shared genetic architecture between traits (genetic correlation).

Key publications:

- Bulik-Sullivan, B., Loh, PR., Finucane, H. et al. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nat Genet 47, 291–295 (2015). https://doi.org/10.1038/ng.3211
- Bulik-Sullivan, B., Finucane, H., Anttila, V. et al. An atlas of genetic correlations across human diseases and traits. Nat Genet 47, 1236–1241 (2015). https://doi.org/10.1038/ng.3406
- Finucane, H., Bulik-Sullivan, B., Gusev, A. et al. Partitioning heritability by functional annotation using genome-wide association summary statistics. Nat Genet 47, 1228–1235 (2015). https://doi.org/10.1038/ng.3404

## Installation

### Prerequisites

- **Python**: Version >3.10
- **Git**: For cloning the repository
- **Poetry**: Python dependency management tool

### Steps

1. **Clone the Repository**

   ```bash
   git clone https://github.com/abdenlab/ldsc-python3.git
   cd ldsc-python3
   ```

2. **Install Poetry**

   If you don't have Poetry installed, you can install it using the following command:

   ```bash
   # From python-poetry.org
   curl -sSL https://install.python-poetry.org | python3 -
   # Or install with pip
   pip install poetry
   ```

   Make sure to add Poetry to your PATH as instructed after installation.

3. **Install Dependencies**

   Use Poetry to install all project dependencies:

   ```bash
   poetry install
   ```

   This will create a virtual environment and install all required packages as specified in `pyproject.toml`.

4. **Activate the Virtual Environment**

   ```bash
   poetry shell
   ```

5. **Verify Installation**

   Run the help command to verify that LDSC is installed correctly:

   ```bash
   python ldsc.py -h
   python munge_sumstats.py -h
   ```

   If these commands display help messages with available options, the installation was successful.

## Running LDSC

LDSC provides several functionalities, including estimating LD Scores, heritability, partitioned heritability, genetic correlation, and the LD Score regression intercept.

### Estimating LD Scores

```bash
ldsc --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC. --l2 --ld-wind-cm 1 --out eur_ldscores
```

### Estimating Heritability

```bash
ldsc --h2 sumstats.txt --ref-ld-chr eur_ldscores/ --w-ld-chr eur_weights/ --out h2_results
```

### Estimating Genetic Correlation

```bash
ldsc --rg trait1.sumstats.gz,trait2.sumstats.gz --ref-ld-chr eur_ldscores/ --w-ld-chr eur_weights/ --out rg_results
```

### Partitioned Heritability

```bash
ldsc --h2 sumstats.txt --ref-ld-chr eur_ldscores/ --w-ld-chr eur_weights/ --overlap-annot --frqfile-chr eur_frq/ --out partitioned_h2
```

Replace the file paths with your actual data files. For more detailed tutorials and options, refer to the [wiki](https://github.com/abdenlab/ldsc-python3/wiki).

## Testing

We have included a comprehensive test suite to ensure the correctness of LDSC.

### Running Tests

1. **Activate the Virtual Environment**

   ```bash
   poetry shell
   ```

2. **Run Tests with Nose2**

   ```bash
   nose2
   ```

   This will execute all unit tests located in the `test` directory.

### Continuous Integration

We use GitHub Actions for continuous integration. The workflow is defined in `.github/workflows/python-project.yml` and runs tests across multiple Python versions.

## Contributing

We welcome contributions from the community. To contribute:

1. **Fork the Repository**

2. **Create a Feature Branch**

   ```bash
   git checkout -b feature/new-feature
   ```

3. **Make Changes and Commit**

   ```bash
   git commit -am "Add new feature"
   ```

4. **Push to Your Fork**

   ```bash
   git push origin feature/new-feature
   ```

5. **Create a Pull Request**

Please ensure that your code passes all tests and adheres to the project's coding standards before submitting a pull request.

### Coding Standards

- **Formatting**: We use `black` for code formatting.
- **Linting**: Code should pass `flake8` checks.
- **Type Checking**: We use `mypy` for static type checking.
- **Imports**: Organize imports using `isort`.

### Setting Up Development Environment

Install development dependencies:

```bash
poetry install --with dev
```

## Docker Setup

We provide a `Dockerfile` to containerize the application.

### Building the Docker Image

```bash
docker build -t ldsc:dev .
```

### Running the Docker Container

```bash
docker run -it ldsc:dev ldsc -h
```

Adjust the command according to your needs, especially if you have specific scripts to run.

## Citing LDSC

If you use LDSC in your research, please cite the following publications:

- **LD Score Regression Methodology**:

  Bulik-Sullivan, B.K. et al. (2015). *LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies*. Nature Genetics, 47(3), 291–295. [doi:10.1038/ng.3211](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3211.html)

- **Genetic Correlation**:

  Bulik-Sullivan, B. et al. (2015). *An Atlas of Genetic Correlations across Human Diseases and Traits*. Nature Genetics, 47(11), 1236–1241. [doi:10.1038/ng.3406](https://www.nature.com/articles/ng.3406)

- **Partitioned Heritability**:

  Finucane, H.K. et al. (2015). *Partitioning Heritability by Functional Annotation Using Genome-Wide Association Summary Statistics*. Nature Genetics, 47(11), 1228–1235. [doi:10.1038/ng.3404](https://www.nature.com/articles/ng.3404)

- **Continuous Annotation Stratified Heritability**:

  Gazal, S. et al. (2017). *Linkage Disequilibrium–Dependent Architecture of Human Complex Traits Shows Action of Negative Selection*. Nature Genetics, 49(10), 1421–1427. [doi:10.1038/ng.3954](https://www.nature.com/articles/ng.3954)

- **Relation to Haseman-Elston Regression**:

  Bulik-Sullivan, B.K. (2015). *Relationship Between LD Score and Haseman-Elston Regression*. bioRxiv. [doi:10.1101/018283](http://dx.doi.org/10.1101/018283)

## License

This project is licensed under the **GNU General Public License v3.0**. You may obtain a copy of the License at [https://www.gnu.org/licenses/gpl-3.0.en.html](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Authors

- **Brendan Bulik-Sullivan**

  Broad Institute of MIT and Harvard

- **Hilary Finucane**

  MIT Department of Mathematics

- **Thomas Reimonn**

  UMass Chan Medical School

## Support

If you encounter issues or have questions, please consider the following resources:

1. **Wiki and Tutorials**

   Detailed tutorials are available in the [wiki](https://github.com/abdenlab/ldsc-python3/wiki) section.

2. **Frequently Asked Questions**

   Common issues are addressed in the [FAQ](https://github.com/abdenlab/ldsc-python3/wiki/FAQ).

3. **Contact Us**

   For further assistance, you can reach out via the [GitHub Issues](https://github.com/abdenlab/ldsc-python3/issues) page.

## Updating LDSC

To update to the latest version of LDSC:

1. **Navigate to the Project Directory**

   ```bash
   cd ldsc-python3
   ```

2. **Pull the Latest Changes**

   ```bash
   git pull
   ```

   If your local repository is up to date, you will see:

   ```
   Already up-to-date.
   ```

   Otherwise, `git` will display the files that have been updated.

3. **Update Dependencies**

   If dependencies have changed, update them with:

   ```bash
   poetry install
   ```

## Obtaining LD Scores

You can download precomputed LD Scores suitable for various analyses:

- **European LD Scores**: [Download](https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2)
- **East Asian LD Scores**: [Download](https://data.broadinstitute.org/alkesgroup/LDSCORE/eas_ldscores.tar.bz2)

Partitioned LD Scores for heritability estimation are also available [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/).

---

We hope LDSC proves valuable in your research. Your contributions and feedback are greatly appreciated!
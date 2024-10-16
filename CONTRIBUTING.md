# Contributing to `modsem`

Thank you for considering contributing to `modsem`! We welcome contributions to help improve this package for estimating interaction effects in structural equation modeling (SEM). To ensure a smooth collaboration, please follow the guidelines below.

## Getting Started

### Fork and Clone the Repository

1. Fork the repository on GitHub.
2. Clone your fork to your local machine.

```sh
git clone https://github.com/your-username/modsem.git
cd modsem
```

### Setting up your Development Environment

1. Ensure you have R installed on your machine.
2. Install the package dependencies.

```r
install.packages("devtools")
devtools::install_deps()
```

3. Install the `modsem` package from your local repository.

```r
devtools::install()
```

## Making Changes

### Creating a Branch

1. Always create a new branch for your work.

```sh
git checkout -b your-branch-name
```

### Making Your Changes

1. Make your changes in the codebase.
2. Ensure that your changes are well-documented.
3. Write tests for your changes if applicable.

### Contributing to Vignettes

We also encourage contributions to the vignettes. If you have a new use case or example, feel free to add or alter vignettes to help demonstrate the functionality of `modsem`.

### Running Tests

1. Run the tests to ensure your changes do not break existing functionality.

```r
devtools::test()
```

## Submitting Your Changes

1. Push your changes to your fork.

```sh
git push origin your-branch-name
```

2. Open a pull request on GitHub against the `main` branch of the original repository.

### Pull Request Guidelines

1. Provide a clear and descriptive title for your pull request.
2. Describe the changes you made and why they are necessary.
3. Reference any related issues or pull requests.
4. Ensure all tests pass and there are no merge conflicts.

## Reporting Issues

If you encounter any issues or have suggestions for improvements, please open an issue on GitHub. Provide as much detail as possible to help us understand and address the issue.
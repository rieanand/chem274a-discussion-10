# Discussion 10 Workbook: Why Does CO₂ Warm the Planet but N₂ Doesn't?
### Learning the Eigen linear algebra library

Earth cools itself by radiating infrared (IR) light to space. N₂ and O₂ (99% of the air) let that light pass, but
greenhouse gases absorb it, because the light matches one of their **vibrational frequencies**. Today you'll write a
C++ program that predicts a molecule's vibrational frequencies, and decides whether it absorbs in the
"atmospheric window" (800–1250 cm⁻¹) where Earth's heat normally escapes.

The goal is learning the [Eigen](https://libeigen.gitlab.io/eigen/docs-nightly/) library and
some C++ habits that will save you hours of debugging.

**Setup** (run once):

```bash
curl -LO https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
tar xjf eigen-3.4.0.tar.bz2
g++ -std=c++17 -Wall -I eigen-3.4.0 main.cpp -o main && ./main
```

**0.1** Eigen is "header-only." Looking at the compile command, what does that mean for you?

---

## Exercise 1: Storing molecules in matrices

To simulate a gas molecule, we first need its atoms' positions. `main.cpp` stores water (the biggest natural
greenhouse gas) and formaldehyde (an air pollutant) as `Eigen::MatrixXd`, with one row per atom and columns x, y, z.
Rotating a molecule is a matrix multiply. Add:

```c++
Eigen::Matrix3d R{{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}};   // 90° rotation about z
Eigen::MatrixXd rotated = h2o * R.transpose();
std::cout << rotated << "\n";
```

**1.1** `Matrix3d` vs `MatrixXd`: what do the `3`, the `X`, and the `d` mean? 

**1.2** Try it: `std::cout << sizeof(Eigen::Matrix3d) << " " << sizeof(Eigen::MatrixXd) << "\n";`
Both hold 9 doubles here. Why is the `MatrixXd` object *smaller*? 

**1.3** Try it: `Eigen::MatrixXd bad = R * form;` (`form` is 4 × 3).
Does it compile? Does it run?
Why can't the compiler catch this mistake? 

**1.4** Try it:

```c++
Eigen::MatrixXd a{{1, 2}, {3, 4}};
a = a.transpose();
std::cout << a << "\n";
```

Why does this not work?
Now add `-DNDEBUG` to the compile command (that's how optimized "release" builds are made) and fix the error. 

---

## Exercise 2: Choosing an eigensolver

A molecule's vibrational frequencies come from the **eigenvalues** of a matrix. Before we trust Eigen with real
data, let's test its two solvers on a practice matrix, `m_to_diag` in `main.cpp`:

```c++
Eigen::EigenSolver<Eigen::MatrixXd> general(m_to_diag);
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> symm(m_to_diag);
std::cout << general.eigenvalues() << "\n\n" << symm.eigenvalues() << "\n";
```

**2.1** What is the `<Eigen::MatrixXd>` part called in C++, and why does the solver need it? 

**2.2** The two solvers give different answers. Look up
[`SelfAdjointEigenSolver`](https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1SelfAdjointEigenSolver.html).
Which part of the matrix does it actually read? Does it check that your matrix is symmetric? 

**2.3** Write a symmetry check and run it:

```c++
std::cout << (m_to_diag - m_to_diag.__________()).cwiseAbs().maxCoeff() << "\n";
```

Which solver's answer is right? 

**2.4** Try it: `Eigen::VectorXd ev = general.eigenvalues();`
What type should `ev` be? 

**2.5** From the docs: are `general`'s eigenvalues sorted? Are `symm`'s? 

---

## Exercise 3: Reading data from a file

The matrix that describes how stiff a molecule's bonds are (the "Hessian") comes from quantum chemistry software.
Each file holds one matrix, written row by row, with 3 rows and 3 columns per atom (x, y, z). Fill in the blanks:

```c++
Eigen::MatrixXd read_hessian(std::string path, int natoms)
{
    Eigen::MatrixXd H(____ * natoms, ____ * natoms);
    std::ifstream in(path);
    for (int i = 0; i < H.rows(); i++)
        for (int j = 0; j < H.cols(); j++)
            in >> H(i, j);
    return H;
}
```

```c++
Eigen::MatrixXd H = read_hessian("formaldehyde_hessian.dat", 4);
std::cout << H.rows() << " x " << H.cols() << "\n" << H << "\n";
```

**3.1** Formaldehyde's matrix is `_____ x _____`.

**3.2** `formaldehyde_hessian.dat` has 12 numbers per line, but `water_hessian.dat` has 3. Does the loop care? Why
or why not? 

**3.3** Try it: change `H(i, j)` to `H[i][j]`. The compiler error contains a message in ALL CAPS. What does it tell
you?

---

## Exercise 4: Loops vs. Eigen array operations

Heavy atoms vibrate more slowly than light ones, so each matrix element is divided by the square root of two atoms'
masses. Formaldehyde's atoms are C, O, H, H.

**Loop version:**

```c++
std::vector<double> mass = {atomic_weights.at("C"), atomic_weights.at("O"),
                            atomic_weights.at("H"), atomic_weights.at("H")};
Eigen::MatrixXd Hmw(H.rows(), H.cols());
for (int i = 0; i < H.rows(); i++)
    for (int j = 0; j < H.cols(); j++)
        Hmw(i, j) = H(i, j) / std::sqrt(mass[i / 3] * mass[j / 3]);
```

**4.1** Rows 0–2 are atom 0, rows 3–5 are atom 1, and so on. For `i = 7`, `i / 3` is `_____`. Why does this only work
because `i` is an `int`? 

**4.2** `atomic_weights.at("C")` vs `atomic_weights["C"]`: what's the difference if the key is misspelled?
(Hint: `atomic_weights` is `const`.) 

**No-loop version:**

```c++
Eigen::VectorXd m3(H.rows());
for (int i = 0; i < H.rows(); i++) m3(i) = mass[i / 3];
Eigen::MatrixXd Hmw2 = H.array() / (m3 * m3.transpose()).array().sqrt();
std::cout << (Hmw - Hmw2).norm() << "\n";
```

**4.3** `m3` is 12 × 1. What shape is `m3 * m3.transpose()`? 

**4.4** What does `.array()` change about how `/` and `*` behave? 

**4.5** The printed difference is `_____`. Do the two versions agree? 

**4.6** For two square matrices `A` and `B`: which one is matrix multiplication, `A * B` or `A.array() * B.array()`? 

---

## Exercise 5: Is formaldehyde a greenhouse gas?

Each eigenvalue gives one vibrational frequency. Let's compute them and see if any land in the atmospheric window.

```c++
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hmw);
if (es.info() != Eigen::__________) throw std::runtime_error("diagonalization failed");

for (int k = 0; k < es.eigenvalues().size(); k++) {
    double lam = es.eigenvalues()(k);
    double nu = 5140.484532 * std::sqrt(std::abs(lam));   // unit conversion to cm^-1
    if (lam < 0) nu = -nu;
    std::cout << nu << "\n";
}
```

**5.1** Fill in the `info()` blank (check the docs). 

**5.2** Use your check from 2.3 on `Hmw`. Result: `_____`. So is `SelfAdjointEigenSolver` safe to use here? 

**5.3** Try it: remove `std::abs`. What prints for the negative eigenvalues? 

**5.4** Why `es.eigenvalues()(k)` and not `es.eigenvalues()[k]`? Would `[k]` work here? (Compare with 3.3.) 

**5.5** `es.eigenvectors().col(k)` is the eigenvector for eigenvalue `k`. Print the `.norm()` of one. You get `_____`.

**Check your output:** ignore the 6 values near zero. Your 6 largest should each be within about 4% of the measured
values 1167, 1249, 1500, 1746, 2783 and 2843 cm⁻¹
([NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C50000&Mask=800)). Pass? 

**5.6** Do any of your frequencies fall in the 800–1250 cm⁻¹ window? 

---

**Bonus:** time the loop version vs. the `.array()` version of Exercise 4 on a 3000 × 3000 random matrix
(`Eigen::MatrixXd::Random(3000, 3000)`) using `<chrono>`, with and without `-O2`. Which is faster? Why?

# *De Novo* Genome Assembly

## The Greedy Shortest Common Superstring
### Task
* In R, implement a function `GreedySuperstring()` according to the pseudocode.

* Input:
    * `S` A `DNAStringSet` object of strings (reads).

* Output:
   * `S` A `DNAStringSet` object of the shortest common superstring (contig).

> **Hint:** 
> Create also functions:
> * `Overlap()` to calculate overlap between two sequences.
> * `OverlapMatrix()` to create a matrix of overlaps among all sequences in `S`.

```
GreedySuperstring(S)
1   while length of S > 1
2     overlapMat <- OverlapMatrix(S)
3     if max(overlapMat) = 0
4       return S
5     else
6       seq1, seq2 ← Two sequences from S with the longest overlap
7       Merge seq1 and seq2 and add the new sequence to S
8       Remove seq1 and seq2 from S
9   return S
```

library(Biostrings)

#' Vypočítá délku nejdelšího překryvu (overlap) mezi sekvencemi.
#' @param seq1 První DNA sekvence (Biostrings::DNAString).
#' @param seq2 Druhá DNA sekvence (Biostrings::DNAString).
#' @return Délka nejdelšího překryvu (integer).
Overlap <- function(seq1, seq2) {
  len1 <- length(seq1)
  len2 <- length(seq2)
  max_overlap <- 0

  # Procházíme všechny možné délky překryvu, od nejdelší po nejkratší (min. 1)
  for (k in min(len1, len2):1) {
    # Zkontroluje, zda konec seq1 o délce k se shoduje se začátkem seq2 o délce k
    if (subseq(seq1, start = len1 - k + 1, end = len1) == subseq(seq2, start = 1, end = k)) {
      max_overlap <- k
      break # Nalezli jsme nejdelší překryv, můžeme skončit
    }
  }
  return(max_overlap)
}

#' Vytvoří matici překryvů pro sadu DNA sekvencí.
#' @param S DNAStringSet objekt sekvencí (čtení).
#' @return Matice (matrix) délek překryvů.
OverlapMatrix <- function(S) {
  n <- length(S)
  overlapMat <- matrix(0, nrow = n, ncol = n)
  rownames(overlapMat) <- names(S)
  colnames(overlapMat) <- names(S)

  # Procházíme všechny dvojice (i, j) kde i != j
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        overlapMat[i, j] <- Overlap(S[[i]], S[[j]])
      }
    }
  }
  return(overlapMat)
}

#' Sestaví Superstring pomocí hladového algoritmu (Greedy Shortest Common Superstring).
#' @param S DNAStringSet objekt sekvencí (čtení).
#' @return DNAStringSet objekt obsahující výsledný kontig (nebo kontigy).
GreedySuperstring <- function(S) {
  if (length(S) == 0) {
    return(DNAStringSet())
  }
  
  # Přiřazení názvů, pokud chybí
  if (is.null(names(S))) {
    names(S) <- paste0("seq", 1:length(S))
  }
  
  # Hladová smyčka
  while (length(S) > 1) {
    # Krok 2: Vytvoř OverlapMatrix
    overlapMat <- OverlapMatrix(S)
    
    # Krok 3: Najdi maximální překryv
    max_overlap <- max(overlapMat)
    
    # Krok 3-4: Pokud je max překryv 0, nelze už dále spojovat
    if (max_overlap == 0) {
      # Vrátí sadu nespojených sekvencí
      return(S) 
    } else {
      # Krok 6: Najdi indexy sekvencí s nejdelším překryvem
      max_indices <- which(overlapMat == max_overlap, arr.ind = TRUE)
      
      # Vezmi první dvojici s největším překryvem (indexy v rámci matice)
      row_idx <- max_indices[1, 1]
      col_idx <- max_indices[1, 2]
      
      # Získání sekvencí (pomocí jmen pro bezpečnější manipulaci)
      seq1_name <- rownames(overlapMat)[row_idx]
      seq2_name <- colnames(overlapMat)[col_idx]
      
      seq1 <- S[[seq1_name]]
      seq2 <- S[[seq2_name]]
      
      # Krok 7: Sloučení sekvencí
      # Sloučená sekvence je seq1 + zbytek (non-overlap část) seq2
      # Délka non-overlap části seq2 je délka(seq2) - max_overlap
      merged_seq <- xscat(seq1, subseq(seq2, start = max_overlap + 1))
      
      # Nové jméno kontigu
      new_name <- paste0(seq1_name, "_", seq2_name, "_merged")
      
      # Přidání nové sekvence do S
      S[[new_name]] <- merged_seq
      
      # Krok 8: Odstranění původních sekvencí
      S <- S[names(S) != seq1_name & names(S) != seq2_name]
    }
  }
  
  # Krok 9: Vrať výsledný (Super)string
  return(S)
}

# Příklad vstupních DNA sekvencí (čtení)
reads_data <- c(
  "ATGCGTAGCT", # seq1
  "TAGCTAGCAT", # seq2 - překrývá se s seq1 (TAGCT)
  "GCATCGATT"   # seq3 - překrývá se s seq2 (GCAT)
)
names(reads_data) <- paste0("read", 1:3)

# Vytvoření DNAStringSet objektu
S_input <- DNAStringSet(reads_data)
print(S_input)

# Testování OverlapMatrix
mat <- OverlapMatrix(S_input)
print(mat)
# Očekávané výsledky:
# read1 -> read2: overlap je 5 (TAGCT)
# read2 -> read3: overlap je 4 (GCAT)

# Spuštění hlavního algoritmu
final_contigs <- GreedySuperstring(S_input)

cat("\n--- Výsledný Superstring ---\n")
print(final_contigs)
# Očekávaný výsledek: ATGCGTAGCATCGATT

<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Basic Git settings</summary>

>* Configure the Git editor
>    ```bash
>    git config --global core.editor notepad
>    ```
>* Configure your name and email address
>    ```bash
>    git config --global user.name "Zuzana Nova"
>    git config --global user.email z.nova@vut.cz
>    ```
>* Check current settings
>    ```bash
>    git config --global --list
>    ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Clone forked repository from your GitHub page to your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_09.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create a new commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send a new commit to your GitHub repository.
```bash
git push origin main
```

</details>

# **MIDS (Multi-frame Individual Distinctions by Size) version 1.0**

**For troubleshooting contact:** [Hayden.Swift@newcastle.edu.au](mailto:Hayden.Swift@newcastle.edu.au)

---

## **MIDS guide on how to implement the analysis and software package in RStudio.**

MIDS package provides tools for applying **Multi-frame Individual Distinctions by Size** to stereo-video camera data.

The package estimates the number of **distinct individuals** observed within a camera deployment using **stasitical inference** and by applying the **MIDS logic** to remove repeat observations across time.

---

## **Package Installation**

```r
install.packages("remotes")   # if not already installed

remotes::install_github("HSUON/MIDS")
```

After installation, MIDS can be loaded with:

```r
library(MIDS)
```

### **MIDS package functions**

**Function 1: `run_mids`** - This function is the main MIDS analysis where individuals are statisically compared using the MIDS logic.

**Function 2: `summarise_mids`** - This function provides the total number of observations from the original dataset (pre-mids), the total number from the MIDS filtered dataset, and the total number of removal individuals.

**Function 3: `summarise_mids_by_group`** - This function provides a detailed list of the individuals that were both retained and removed from the analysis.

---

## **Required columns for the analysis:**

### **MIDS requires:**

- a grouping structure, which includes a **"Camera"** number/ID coloumn and a **"Species"** coloumn. However, other factors such as **"location"**, **"season"** or **"site"** can be added if required to group the cameras.
- a **"frame"** or time column for when the individual(s) has/have appeared within a video file
- a **"length"** column with the size of each individual (mm)
- a **"precision_error"** column including the measurement error for each individual size measurement (mm)

---

## **Example test data**

```r
test_data <- data.frame(
Camera = c("C1","C1","C1","C1","C2","C2"),
Frame = c(1,1,2,2,1,2),
Species = c("A","A","A","A","A","A"),
Length = c(105,112,109,150,98,104),
Precision_Error = c(2.5,2.0,3.0,2.5,1.8,2.2)
)
```

## **Example of how to write and run the MIDS functions:**

### **Main MIDS analysis (ensure any other grouping factors are included within "group_cols")**

```r
results <- run_mids(
data = test_data,
group_cols = c("Camera", "Species"),
frame_col = "Frame",
length_col = "Length",
precision_col = "Precision_Error"
)

results
```

## **Summary of total observations**

```r
summarise_mids(test_data, results)
```

## **Summary of individuals retained and removed by group (ensure any other grouping factors are included within "group_cols")**

```r
summarise_mids_by_group(
original_data = test_data,
mids_results = results,
group_cols = c("Camera", "Species")
)
```

## License

MIDS is made available under the **PolyForm Noncommercial License 1.0.0**.

The software may be used, modified, and distributed for permitted non-commercial purposes in accordance with the license terms. This includes use by educational institutions, public research organisations, environmental protection organisations, and government institutions.

**Commercial use of MIDS is not permitted under this license and requires a separate commercial license from the copyright holder.**

For commercial licensing enquiries, please contact:  
**Hayden Swift** — Hayden.Swift@newcastle.edu.au

See the [LICENSE](LICENSE) file for the complete license terms.

© 2026 Hayden Swift.

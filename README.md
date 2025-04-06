
# 🧬 Protein Structure Prediction using Genetic Algorithm (MATLAB GUI)

This project demonstrates a **simple yet interactive Genetic Algorithm (GA)** for predicting **2D protein folding structures**. Built in **MATLAB**, the GUI allows users to tune GA parameters and visually track the evolution of protein folding toward an optimal configuration.

🎥 _A video demonstration is available — check it out to see the protein structure evolve live!_

---

## 🚀 Features

- ✅ **User-friendly MATLAB GUI**
- 🎛️ Adjustable GA parameters:
  - Population size
  - Number of genes (protein chain length)
  - Mutation rate
  - Crossover rate
  - Number of generations
- 📈 Live plotting of **best fitness over generations**
- 🧩 Real-time 2D visualization of the **protein structure evolving**
- ⏳ Animated progression with delay to help visualize how the structure improves

---

## 📌 How It Works

Each individual in the population is a vector of angles, representing fold directions of the protein chain. The goal is to **minimize a fitness function** (e.g., sum of squared angles) to find a stable configuration.

### Visualization
- The folding starts from the origin and builds the chain step by step using angles.
- You’ll see the folding structure update **in each generation**, reflecting the best configuration so far.

---

## 🧠 Fitness Function

```matlab
fVal = sum(individual.^2);
```

A simple sum-of-squares is used for illustration. This can be extended to more biologically accurate models (e.g., hydrophobic-hydrophilic models, energy-based models).

---

## ▶️ How to Run

1. Open `ProteinFoldingGA_GUI.m` in MATLAB.
2. Click **Run** or type `ProteinFoldingGA_GUI` in the Command Window.
3. Tune the parameters in the GUI.
4. Click **Start GA** to watch the magic unfold!

---

## 📹 Video Demo

> 🎬 https://youtu.be/OXN7CNdOnkw 

---

## 🧪 Future Improvements

- 3D protein folding visualization
- More realistic energy-based fitness functions
- Real protein sequence input support
- Export final structure to PDB format

---

## 👨‍💻 Created By

**Roshan Rateria**  
_Exploring the intersection of bioinformatics & computation 🤖🧬_

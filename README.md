# 🌌 CosmoGraph: The Computational Universe Project

**CosmoGraph** est une simulation cosmologique expérimentale basée sur la théorie des graphes.  
Au lieu de simuler des particules dans un espace pré-défini, ce projet génère l'espace-temps lui-même comme un réseau dynamique, démontrant comment les lois physiques (gravité, expansion, flèche du temps) peuvent émerger de la topologie pure.

---

## 🔬 La Théorie

Le projet repose sur l'hypothèse que l'univers est un graphe sans échelle (Scale-Free Network).

1. **Genèse :** Utilisation de l'algorithme Barabási-Albert (attachement préférentiel).  
2. **Matière :** Les nœuds à haute centralité (Hubs) représentent la matière massive.  
3. **Matière Noire :** Les nœuds périphériques à faible degré assurent la cohésion du graphe.  
4. **Énergie Noire :** L'ajout constant de nouveaux nœuds étire les géodésiques (Expansion de l'univers).  
5. **Flèche du Temps :** L'entropie de Von Neumann du réseau augmente de manière monotone ($dS/dt > 0$).
6. **Relativité discrète :** Masse et courbure corrèlent fortement (Forman-Ricci, corr. ≈ -0.898).

---

## 🚀 Fonctionnalités du Dashboard

Interface interactive (Streamlit) permettant d'explorer :

- **Structure 3D** : Visualisation rotative de l'univers et de ses amas.  
- **Gravité & preuve** : Orbite entropique et preuve 1/r².  
- **Matière Noire** : Cartographie des anomalies topologiques (haute centralité / faible degré).  
- **Expansion** : Simulation de la loi de Hubble (Distance moyenne vs Volume).  
- **Multivers** : Scan Small-World (σ) montrant la zone habitable (m≈2).  
- **Entropie** : Flèche du temps (S_VN(t) croissante).  
- **Relativité** : Corrélation Masse/Courbure (Forman-Ricci).  
- **Monde Réel (Toile cosmique)** : Carte de courbure sur une toile cosmique simulée (RGG).

---

## 🛠 Installation et Utilisation

### Prérequis
- Python 3.8+
- Environnement virtuel recommandé (à créer localement, non fourni dans le dépôt)

### Installation

```bash
# Cloner le dépôt
git clone https://github.com/LudovicTagnon/CosmoGraph-Emergent-Gravity.git
cd CosmoGraph-Emergent-Gravity

# Créer un environnement virtuel
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Installer les dépendances
pip install -r requirements.txt
```

### Lancer la Simulation

Pour régénérer les données et lancer le Dashboard interactif :

```bash
# (Optionnel) Recalculer les jeux de données / figures clés
cd python_project
python scripts/multiverse_scan.py
python scripts/arrow_of_time.py
python scripts/detect_dark_matter.py
python scripts/simulate_expansion.py
python scripts/power_spectrum_combined.py
python scripts/fetch_and_test_real_sdss.py
python scripts/sensitivity_analysis.py
python scripts/plot_sensitivity.py

# Lancer l'interface
streamlit run app.py
```

### Données externes
- Les scripts utilisent un échantillon SDSS local. Vous pouvez tester d’autres catalogues publics, par ex. le dataset Kaggle “SDSS Galaxy Classification DR18” : https://www.kaggle.com/datasets/bryancimo/sdss-galaxy-classification-dr18 (non inclus dans le dépôt). Placez le CSV dans `python_project/data/` et adaptez `fetch_and_test_real_sdss.py` si besoin.

---

## 📊 Résultats Clés

| Métrique | Valeur Calculée | Interprétation |
| :--- | :--- | :--- |
| **Coefficient Small-World ($\sigma$)** | 5.64 (à $m=2$) | Structure optimale pour le transfert d'info |
| **Entropie Finale** | 4.30 | Flèche du temps (S(t) croissante) |
| **Loi de Hubble** | Confirmée | Distance/Volume augmente |
| **Einstein Score (Forman-Ricci)** | ≈ -0.898 | Masse et courbure fortement corrélées |
| **Matière Noire** | ~5% d'anomalies | Gravité fantôme via topologie |

---

*Projet généré par une collaboration Homme-IA explorant la cosmologie computationnelle.*  
*Science via la topologie : Gravité, Matière Noire, Expansion, Flèche du temps, Relativité émergent d'un graphe simple.* 

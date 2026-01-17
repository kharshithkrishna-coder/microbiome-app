# Gut Microbiome Analysis Dashboard

A shareable web application for analyzing gut microbiome composition and predicting nutrient bioavailability.

## Features

- 🧬 **OTU Analysis** - Visualize microbial species abundance and distribution
- 📊 **Nutrient Absorption Modeling** - Predict bioavailability of 6 key nutrients
- 🔬 **Interactive Simulations** - Test the impact of microbiome perturbations
- 📈 **Dynamic Visualizations** - Plotly-based interactive charts
- 🌐 **Shareable Web Interface** - Deploy and share results easily

## Installation

### 1. Install Dependencies

```bash
pip install streamlit pandas numpy scipy plotly
```

### 2. Prepare Data

Ensure `india_species_abundance_clean.tsv` is in the same directory as `app.py`.

## Running Locally

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

## Deployment Options

### Option 1: Streamlit Cloud (Recommended)

1. Push your code to GitHub:
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin <your-repo-url>
   git push -u origin main
   ```

2. Go to [share.streamlit.io](https://share.streamlit.io)

3. Sign in with GitHub and deploy:
   - Click "New app"
   - Select repository, branch, and file (`app.py`)
   - Click "Deploy"

Your app will be live at `https://<your-username>-<app-name>.streamlit.app`

### Option 2: Heroku

1. Create `requirements.txt`:
   ```bash
   pip freeze > requirements.txt
   ```

2. Create `Procfile`:
   ```
   web: streamlit run app.py --logger.level=error
   ```

3. Deploy:
   ```bash
   heroku login
   heroku create <app-name>
   git push heroku main
   ```

### Option 3: Docker

1. Create `Dockerfile`:
   ```dockerfile
   FROM python:3.9
   WORKDIR /app
   COPY requirements.txt .
   RUN pip install -r requirements.txt
   COPY . .
   EXPOSE 8501
   CMD ["streamlit", "run", "app.py"]
   ```

2. Build and run:
   ```bash
   docker build -t microbiome-app .
   docker run -p 8501:8501 microbiome-app
   ```

## File Structure

```
project/
├── app.py                              # Main Streamlit app
├── gut_model.ipynb                     # Jupyter notebook version
├── gut_model.py                        # Python script version
├── india_species_abundance_clean.tsv   # OTU abundance data
├── requirements.txt                    # Python dependencies
└── README.md                          # This file
```

## Sections

### 1. Overview
- High-level statistics and key metrics
- Top 10 most abundant species visualization

### 2. OTU Analysis
- Species abundance distribution
- Cumulative abundance curves
- Complete species data table

### 3. Nutrient Absorption
- Baseline absorption scores for all nutrients
- Nutrient-specific trait contributions
- Interactive nutrient selection

### 4. Simulations
- Perturb individual bacterial species
- Simulate abundance changes
- Real-time impact on nutrient absorption

### 5. About
- Project methodology
- Data requirements
- Technologies used

## Data Format

Input OTU table should be formatted as:
- Rows: Species names (e.g., `Lactobacillus_plantarum`)
- Columns: Sample IDs
- Values: Relative abundance counts

Example:
```
                    Sample1  Sample2  Sample3
Lactobacillus_plantarum   1000     1200     900
Bifidobacterium_longum     800     1100     950
Bacteroides_vulgatus      1500     1400     1300
```

## Nutrient Model

The model includes 6 key nutrients:
- **Iron** - Affected by SCFA, pH reduction, barrier support, siderophore production
- **Vitamin B12** - Requires vitamin biosynthesis, barrier support
- **Folate** - Primary from vitamin biosynthesis
- **Calcium** - Enhanced by SCFA production, pH reduction
- **Magnesium** - Driven by SCFA production
- **Zinc** - Requires barrier support, hindered by siderophore competition

## Functional Traits

Six key bacterial traits are modeled:
1. **SCFA Production** - Short-chain fatty acid fermentation
2. **pH Reduction** - Acidification of intestinal environment
3. **Barrier Support** - Tight junction maintenance
4. **Vitamin Biosynthesis** - Endogenous vitamin production
5. **Siderophore Production** - Iron sequestration capability

## Sharing Tips

- **Full URL**: Share the deployed app URL directly
- **Screenshots**: Use the "Share" button on Streamlit Cloud
- **Embed**: Embed the app in a website using an iframe
- **Export Data**: Download simulation results from the dashboard

## Troubleshooting

**"FileNotFoundError: india_species_abundance_clean.tsv"**
- Ensure the data file is in the same directory as `app.py`
- Check file name spelling exactly

**"ModuleNotFoundError: No module named 'streamlit'"**
- Run: `pip install streamlit`

**App runs slowly**
- Increase `@st.cache_data` timeout
- Consider caching more functions
- Run simulations with fewer bootstraps

## Future Enhancements

- [ ] User file upload support
- [ ] Advanced network visualization (Cytoscape)
- [ ] Machine learning predictions
- [ ] Publication-ready export formats
- [ ] Multi-user support with database backend

## License

This project is part of an academic research initiative.

## Contact

For questions or suggestions, please reach out to the research team.

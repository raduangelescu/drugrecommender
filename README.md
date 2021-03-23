# GeneDrugRecommender
- Recommends drug combinations for a certain bad gene expression
- Tries to figure out what disease it is based on a certain bad gene expression

# How to use
- Have Python 3.9 and install **requirements.txt** using pip
- Run startup.py
- Quick start
```python
from genedrugrecommender.drugrecommender import DrugRecommender

recommender = DrugRecommender('data/config.json')
recommender.run('GSE15852')
```
- You may change **config.json** 

from genedrugrecommender.drugrecommender import DrugRecommender

recommender = DrugRecommender('data/config.json')
recommender.run('GSE15852')

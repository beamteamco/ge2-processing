from sklearn.cluster import DBSCAN
import numpy as np
#data = np.random.rand(500,3)
data = np.array([[1, 2, 1],
		[1, 2, 1.1],
		[0.9, 2, 1]])

db = DBSCAN(eps=0.12, min_samples=1).fit(data)
labels = db.labels_
from collections import Counter
print Counter(labels)


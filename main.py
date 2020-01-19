import numpy as np
from Bio import SeqIO
from bioenhancers.features import KmerVectorizer
from sklearn.ensemble import RandomForestClassifier

vectorizer: KmerVectorizer = KmerVectorizer(k=4, step=1)
positive_seqs = [str(record.seq).upper()
                 for record in SeqIO.parse("./data/vista1500", "fasta")]
negative_seqs = [str(record.seq).upper()
                 for record in SeqIO.parse("./data/randoms1500", "fasta")]


positive_features = vectorizer.transform(positive_seqs)
negative_features = vectorizer.transform(negative_seqs)
positive_features_labs = np.ones(positive_features.shape[0], )
negative_features_labs = np.zeros(negative_features.shape[0], )

features = np.vstack((positive_features, negative_features))
labs = np.hstack((positive_features_labs, negative_features_labs))

rf_classifier = RandomForestClassifier()
rf_classifier.fit(X=features, y=labs)

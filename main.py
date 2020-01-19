import pdb
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

chr21 = SeqIO.read("./data/chr21.fa", "fasta")
chr21_seq = str(chr21.seq).upper()
chunk_size = 1500
chunk_step = 750
chr21_chunks = [chr21_seq[ind:(ind+chunk_size)]
                for ind in range(0, len(chr21_seq), chunk_step)]
chr21_n_chunks_ind = [ind for ind,
                      chunk in enumerate(chr21_chunks) if 'N' in chunk]
chr21_clean_chunks_ind = [ind for ind,
                          chunk in enumerate(chr21_chunks) if 'N' not in chunk]

chr21_features = vectorizer.transform(
    [chr21_chunks[ind] for ind in chr21_clean_chunks_ind])
enhancer_probas = [enhancer_proba for _,
                   enhancer_proba in rf_classifier.predict_proba(X=chr21_features)]
chr21_probas = [0 for ind in range(len(chr21_chunks))]
mean_proba = np.mean(enhancer_probas)

chr21_probas_all = np.zeros((len(chr21_chunks, )))
chr21_probas_all[np.array(chr21_n_chunks_ind)] = mean_proba
chr21_probas_all[np.array(chr21_clean_chunks_ind)] = enhancer_probas

with open("data/chr21.wig", "w") as f:
    f.write(
        f"fixedStep chrom=chr21 start=0 step={chunk_step} span={chunk_size}\n")
    f.write("\n".join([str(x) for x in chr21_probas_all.tolist()]))

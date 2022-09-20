.. _association:

===========================
Association
===========================


InfraPy uses a pair-based, joint-likelihood method to identify those detection sets with high likelihood of originating from a common source.  Events are identified using a pair-based Bayesian algorithm that defines the association between detection pairs from their joint-likelihood and identifies events via hierarchical clustering analysis

- The current implementation utilizes only spatial and temporal coincidence, but additional detection information can further improve event identification

- Evaluation using a synthetic data set has shown some mixing of spatially similar events poorly resolved by network geometry and occasional inclusion of "noise" detections in event clustering

See the following for more references on Association:

- `Blom et al., 2020 <https://academic.oup.com/gji/advance-article-abstract/doi/10.1093/gji/ggaa105/5800992>`_

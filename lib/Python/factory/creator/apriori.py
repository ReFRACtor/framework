import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class CreatorApriori(Creator):

    apriori = param.Array(dims=1)
    covariance = param.Array(dims=1, required=False)
    initial_guess = param.Array(dims=1, required=False)
    retrieved = param.Scalar(bool, required=False)

    def initial_guess_or_apriori(self):
        ig_val = self.initial_guess()

        if ig_val is None:
            return self.apriori()
        else:
            return ig_val

    def retrieval_flag(self):
        ap = self.apriori()
        retrieved = self.retrieved()

        if retrieved is None or retrieved:
            return np.ones(ap.shape[0], dtype=bool)
        else:
            return np.zeros(ap.shape[0], dtype=bool)

    def sv_initial_guess(self):

        ig = self.initial_guess_or_apriori()
        ap = self.apriori()
        cov = self.covariance()
        rflag = self.retrieval_flag()

        sv_ig = rf.InitialGuessValue()
        sv_ig.apriori_subset(rflag, ap)
        sv_ig.apriori_covariance_subset(rflag, cov)
        sv_ig.initial_guess_subset(rflag, ig)
        
        return sv_ig

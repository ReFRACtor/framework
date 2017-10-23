import numpy as np

from .base import Creator
from .. import param

from refractor import framework as rf

class CreatorApriori(Creator):

    apriori = param.Choice(param.Array(dims=1), param.ArrayWithUnit(dims=1))
    covariance = param.Array(dims=1, required=False)
    initial_guess = param.Array(dims=1, required=False)
    retrieved = param.Scalar(bool, required=False)

    def initial_guess_or_apriori(self):
        ig_val = self.initial_guess()

        if ig_val is None:
            ap = self.apriori()
            if hasattr(ap, "value"):
                return ap.value
            else:
                return ap
        else:
            return ig_val

    def retrieval_flag(self):
        ap = self.apriori()
        retrieved = self.retrieved()

        if hasattr(ap, "value"):
            ap_shape = ap.value.shape
        else:
            ap_shape = ap.shape

        if retrieved is None or retrieved:
            return np.ones(ap_shape, dtype=bool)
        else:
            return np.zeros(ap_shape, dtype=bool)

    def sv_initial_guess(self):

        ig = self.initial_guess_or_apriori()
        ap = self.apriori()
        cov = self.covariance()
        rflag = self.retrieval_flag()

        if hasattr(ap, "value"):
            ap_val = ap.value
        else:
            ap_val = ap

        sv_ig = rf.InitialGuessValue()
        sv_ig.apriori_subset(rflag, ap_val)
        sv_ig.apriori_covariance_subset(rflag, cov)
        sv_ig.initial_guess_subset(rflag, ig)
        
        return sv_ig

class CreatorAprioriMultiChannel(Creator):

    apriori = param.Choice(param.Array(dims=2), param.ArrayWithUnit(dims=2))
    covariance = param.Array(dims=2, required=False)
    initial_guess = param.Array(dims=2, required=False)
    retrieved = param.Iterable(required=False)

    def initial_guess_or_apriori(self):
        ig_val = self.initial_guess()

        if ig_val is None:
            ap = self.apriori()
            if hasattr(ap, "value"):
                return ap.value
            else:
                return ap
        else:
            return ig_val

    def retrieval_flag(self):
        ap = self.apriori()
        retrieved = self.retrieved()

        if hasattr(ap, "value"):
            ap_shape = ap.value.shape
        else:
            ap_shape = ap.shape

        if retrieved is None:
            return np.ones(ap_shape, dtype=bool)
        else:
            flags = np.empty(ap_shape, dtype=bool)
            for chan_idx, chan_flag in enumerate(retrieved):
                flags[chan_idx, :] = chan_flag
            return flags

    def sv_initial_guess(self):
        raise NotImplementedError("Not yet done")

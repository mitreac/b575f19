from collections import Counter, namedtuple

Consensus = namedtuple('Consensus', ('base', 'count'))

# Pileup Object
class Pileup:
    """An object that represents the observed bases and their counts at a specific position
    
    Attributes:
        depth_offset (int): how much depth should be offset for normalization (default: 0)
        counts (collections.Counter): Observed bases and their number of occurrences (default: None)
        depth (int): Sum of all observed base counts
        consensus (collections.namedtuple or None): The most common base and its number of occurrences
        maf (float or None): the mean allele frequency of the consensus base versus depth
    """
    
    __slots__ = 'depth_offset _counts'.split()
    
    def __init__(self, depth_offset = 0, counts = None):
        self.depth_offset = depth_offset
        if counts is None:
            self._counts = Counter()
        else:
            if isinstance(counts, Counter):
                self._counts = counts
            else:
                raise ValueError('counts must be of type collections.Counter')
    
    @property
    def depth(self):
        return sum(self.counts.values()) + self.depth_offset
    
    @property
    def consensus(self):
        """Read-only"""
        consensus = self.counts.most_common(1)
        return Consensus(*consensus[0]) if consensus else None
    
    @property
    def counts(self):
        """Read-only"""
        return self._counts

    @property
    def maf(self):
        """Read-only"""
        if self.consensus:
            return self.consensus.count / self.depth
        else:
            return None

    def update(self, observation):
        """Update the observed base counts
        
        Args:
            observation (str): the observed base at the given position
        
        Usage:
            >>> pile = Pileup()
            >>> pile.update('A')
        """
        self._counts.update(observation)
    
    def __str__(self):
        return f'Pileup(depth = {self.depth + self.depth_offset}, counts = {self.counts}, consensus = {self.consensus}, maf = {self.maf})'
    
    def __repr__(self):
        return f'Pileup({"depth_offset = " + str(self.depth_offset) + ", " if self.depth_offset else ""}counts = {self.counts})'
    
    def __eq__(self, other):
        if isinstance(other, Pileup):
            return self.consensus.base == other.consensus.base
        else:
            raise ValueError('Comparison must be between two Pileup objects')
    
    def __len__(self):
        return len(self.counts)
    
    @classmethod
    def __dir__(cls):
        return [stuff for stuff in dir(cls) if not stuff.startswith('_')]
    
    def to_dict(self):
        return {
            'depth': self.depth,
            'counts': self.counts,
            'consensus': tuple(self.consensus)
        }
    
    def is_valid(self):
        """Check whether or not Pileup object contains data
        
        Returns:
            (bool): True if data present, otherwise False
        """
        return True if self.counts else False
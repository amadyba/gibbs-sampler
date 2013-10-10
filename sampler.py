import random
import sys
LETTERS = ['A', 'C', 'G', 'T']
EPS = 10 ** -10
#TODO: Get letters used in DNA so arbitrary letters can be used.
#TODO: Multiple conversion criteria
#TODO: bug: Every sample gets the same answer!
class MathsIsHardException(Exception):
    pass

class Sampler():

    def __init__(self, filename):

        if len(sys.argv) != 3:
            print """ 
            Usage: (python|pypy) sampler.py <length of motif to search for> <number of identical results needed to acquire "convergence">
            """
            sys.exit(1)

        #we take the length of the motif to find as a command line argument
        self.l = int(sys.argv[1])
        self.k = int(sys.argv[2])
        self.DNA = self.read_DNA(filename)
        # make sure all the sequences are the same length
        assert len(set([len(s) for s in self.DNA])) == 1

        self.t = len(self.DNA)  # number of sequences
        self.n = len(self.DNA[0])  # length of the sequences

        #The list of samples that the sampler takes. One is added each iteration that the sampler hasn't converged.
        self.samples = []

    def check_convergence(self):
        #if the last k samples are identical, we've converged
        last_k_samples = self.samples[-self.k:]
        if len(last_k_samples) == self.k and all(sample == last_k_samples[0] for sample in last_k_samples):
            return True
        return False

    def sample(self):
        starting_positions = self.get_starting_positions()
        self.samples = []
        while not self.check_convergence():
            # generate the lmers that come from the starting positions using quite a few square brackets
            tuples = [self.DNA[i][starting_positions[i]: starting_positions[i] + self.l]
                      for i in xrange(self.t)]
            # print_matrix(DNA)
            # print "Starting sequences: ", tuples

            # choose a sequence from the DNA sequences randomly
            sequence = random.choice(self.DNA)
            sequence_index = self.DNA.index(sequence)
            self.DNA.remove(sequence)
            profile = self.create_profile(tuples)
            #print 'Profile: ', profile

            # for each position i in the chosen DNA sequence, find the
            # probability that the lmer starting in this position is generated
            # by the profile
            probs = []
            for i in xrange(self.n - self.l):
                lmer = sequence[i: i + self.l]
                prob = self.get_generation_prob(lmer, profile)
                probs.append(prob)
            probs = self.normalise(probs)
            #print "probs, ", probs

            new_starting_index = self.choose_from_distribution(probs)
            starting_positions[sequence_index] = new_starting_index

            # don't forget to put the sequence back
            self.DNA.insert(sequence_index, sequence)
            motif = self.get_motif(profile)
            score = self.get_score(profile) * 100 / self.l
            print (motif, score)
            self.samples.append((motif, score))

        #print "Samples: ",self.samples
    
        best_motif = self.max_score(self.samples)
        # print "Best motif:", best_motif
        return best_motif

    def sample_random_starting_points(self, nsamples):
        motifs = []
        for i in xrange(nsamples):
            sample = self.sample()
            motifs.append(sample)
        print "motifs: ,", motifs
        print "Final answer: ", self.max_score(motifs)
    def read_DNA(self, filename):
        return [line.strip() for line in open(filename, 'rU').readlines()]

    def print_matrix(self, dnalist):
        print "\n".join(dnalist)

    def get_starting_positions(self):
        # we can't choose starting positions that cause subsequences that would
        # overflow the length of the sequence
        return [random.randint(0, self.n - self.l) for i in xrange(self.t)]

    def create_profile(self, matrix):
        # letter -> list of normalised probabilities of letter appearing at
        # index i in the matrix
        # 
        # initially maps each nucleotide - > [0, 0, ..... 0] 
        profile = {
            nucleotide: [0 for i in xrange(self.l)] for nucleotide in LETTERS
        }

        #count the occurences of each letter
        for i in xrange(self.t):
            for j in xrange(self.l):
                profile[matrix[i][j]][j] += 1
        # normalise
        for letter, scores in profile.iteritems():
            profile[letter] = [float(i) / self.t for i in scores]
        return profile

    # EVERYBODY STAND BACK THIS PYTHON IS GETTING A BIT INTENSE
    # (This is a bit like a static method in other languages)
    @classmethod
    def normalise(self, dist):
        summation = sum(dist)
        return [x / summation for x in dist]

    def get_generation_prob(self, lmer, profile):
        # the probability of the whole string is the product of the probability
        # of each letter being generated in its position over all letters
        prob = 1  # start at 1 so we can multiply
        for i in xrange(len(lmer)):
            profile_prob = profile[lmer[i]][i]
            # this chance can be 0. We don't want to make the whole total 0, so
            # use a very small value to approximate 0. TODO: is this a hack?
            if profile_prob == 0:
                prob *= EPS
            else:
                prob *= profile_prob
        return prob

    @staticmethod
    def max_score(tup):
        # maximises over an iterable of tuples, in which the 
        # second element in the tuple is the score which to maximise over
        return max(tup, key=lambda x: x[1])

    @staticmethod
    def choose_from_distribution(dist):
        rand = random.random()
        for prob in dist:
            rand -= prob
            if rand <= 0:
                return dist.index(prob)
        raise MathsIsHardException, "Choosing from a distribution is hard. Did you make sure the distribution was normalised?"

    def get_motif(self, profile):
        motif = ""
        for i in xrange(self.l):
            candidates = [(letter, profile[letter][i]) for letter in LETTERS]
            motif += self.max_score(candidates)[0]
        return motif

    # calculate the score of a profile, as defined in lectures
    def get_score(self, profile):
        # maybe these list comprehensions are getting a bit out of hand....
        return sum(max(index_frequencies) for index_frequencies in [[freq[i] for freq in profile.itervalues()] for i in xrange(self.l)])

if __name__ == "__main__":
    sampler = Sampler("test/DNA2.txt")
    sampler.sample_random_starting_points(nsamples=10)

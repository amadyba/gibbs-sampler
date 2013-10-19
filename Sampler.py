import random

import sys
VERY_SMALL = 10 ** -10
# TODO: Multiple conversion criteria


#things the program will say it's doing while it's loading
LOADING_PHRASES = ["Actually finding motifs", 
                   "Reticulating splines", 
                   "Never giving you up",
                   "Never letting you down",
                   "Approximating pi",
                   "Hashing hashmaps",
                   "Moving bytes",
                   "Writing report",
                   "Calculating marking weighting of this program",
                   "Becoming SkyNet",
                   "Studying AI",
                   "Renaming COMP3456 -> BINF3456",
                   "Renaming 'Computational Methods for Life Sciences' -> 'Bioinformatics'",
                   "Recursing",
                   "Cursing",
                   "Re-recursing",
                   "Shining shoes",
                   "Flim-flamming"
                  ]

class MathsIsHardException(Exception):

    """An exception thrown when maths is hard."""
    pass


class Sampler():

    """A class to perform the file I/O and perform the actual gibbs sampling. """

    def __init__(self, filename, l, k):
        """Takes the path to the DNA matrix as a file, with one DNA sequence per line, each of the same length"""

        #the length of the motif to find
        self.l = l
        # the number of iterations for which nothing changing will constitute
        # "convergence"
        self.k = k
        self.DNA = self.read_DNA(filename)
        # make sure all the sequences are the same length
        assert self.all_the_same([len(s) for s in self.DNA])

        # Now, people using this software might want to use letters other than
        # ACTG, so we'll scan to see what letters they DID use.
        self.letters = list(set(
            letter for sequence in self.DNA for letter in sequence))

        self.t = len(self.DNA)  # number of sequences
        self.n = len(self.DNA[0])  # length of the sequences
        if self.l >= self.n:
            print "l was >= n. That doesn't really make sense. Try making l smaller or using a different input file"
            sys.exit(1)
        # The list of samples that the sampler takes. One is added each iteration that the sampler hasn't converged.
        # We'll max over this later
        self.samples = []
        # We'll keep track of the last k profiles, so we can tell when we've
        # converged
        self.profiles = []
        #maps motif - > list of l starting positions
        self.motif_to_starting_positions = {}

    @staticmethod
    def all_the_same(sequence):
        """Returns true if everything in the sequence has the same value as everything else, false otherwie"""
        return all(item == sequence[0] for item in sequence)

    def check_convergence(self):
        last_k_profiles = self.profiles[-self.k:]
        if len(last_k_profiles) < self.k:
            return False
        return self.all_the_same(last_k_profiles)

    def sample(self):
        starting_positions = self.get_starting_positions()
        # reset the lists keeping track of samples and profiles
        self.samples = []
        self.profiles = []
        while not self.check_convergence():
            # generate the lmers that come from the starting positions using
            # quite a few square brackets
            tuples = [self.DNA[i][starting_positions[i]: starting_positions[i] + self.l]
                      for i in xrange(self.t)]
            # choose a sequence from the DNA sequences randomly
            sequence = random.choice(self.DNA)
            # remember where it was in the DNA matrix
            sequence_index = self.DNA.index(sequence)
            # remove it
            self.DNA.remove(sequence)

            # make a new profile from the starting positions
            # we COULD fiddle around adding and subtracting 1 from positions in the old profile, and that would save us an O(tl) computation.
            # But if we don't it makes it cleaner to keep track of profiles.
            # This is because keeping track of each profile requires an O(tl)
            # copying operation anyway, so the worst case running time of this
            # algorithm is the same.
            profile = self.create_profile(tuples)
            # save this profile for convergence checking

            self.profiles.append(profile)
            # if saving this profile overflowed the profile buffer, get rid of
            # the oldest thing in it.
            if len(self.profiles) > self.k:
                self.profiles.pop(0)

            # for each position i in the chosen DNA sequence, find the
            # probability that the lmer starting in this position is generated
            # by the profile
            probs = []
            for i in xrange(self.n - self.l):
                lmer = sequence[i: i + self.l]
                prob = self.get_generation_prob(lmer, profile)
                probs.append(prob)
            probs = self.normalise(probs)
            # print "probs, ", probs

            motif = self.get_motif(profile)
            score = self.get_score(profile) * 100 / self.l
            # print (motif, score)
            self.samples.append((motif, score))
            #save the starting position for later so we can print it
            #this only takes O(l) time (to hash the string)
            self.motif_to_starting_positions[motif] = starting_positions

            new_starting_index = self.choose_from_distribution(probs)
            starting_positions[sequence_index] = new_starting_index

            # don't forget to put the sequence back
            self.DNA.insert(sequence_index, sequence)

        # there's no need to take the sample we converged to as the sample we obtained
        # of all the samples we found, take the best one. Hopefully it's the
        # one we converged to, but it's not always
        best_motif = self.max_score(self.samples)
        return best_motif

    
    def sample_random_starting_points(self, nsamples, silly=False):
        """Start the sampler from nsamples random starting points and take the best answer after that"""
        motifs = []
        #print "Finding motifs ...."
        for i in xrange(nsamples):
            if silly:
                if random.random() < 0.3 and LOADING_PHRASES:
                    phrase = random.choice(LOADING_PHRASES)
                    LOADING_PHRASES.remove(phrase)
                    print phrase + "..."
            sample = self.sample()
            motifs.append(sample)
        motif = self.max_score(motifs)
        #print
        #print "Final answer:", motif[0]
        #print "Score:", motif[1]
        #print 
        #self.show_motif(motif[0])
        return motif

    def read_DNA(self, filename):
        """Read the DNA from the file, converting everything to uppercase for consistency"""
        return [line.strip().upper() for line in open(filename, 'rU').readlines()]

    def print_matrix(self, dnalist):
        print "\n".join(dnalist)

    def get_starting_positions(self):
        # we can't choose starting positions that cause subsequences that would
        # overflow the length of the sequence
        return [random.randint(0, self.n - self.l) for i in xrange(self.t)]

    def create_profile(self, matrix):
        # letter -> list of normalised probabilities of letter appearing at
        # index i in the set of l-tuples represented by matrix
        #
        # initially maps each nucleotide - > [0, 0, ..... 0]
        profile = {
            nucleotide: [0 for i in xrange(self.l)] for nucleotide in self.letters
        }

        # count the occurences of each letter
        for i in xrange(self.t):
            for j in xrange(self.l):
                profile[matrix[i][j]][j] += 1

        # normalise
        for letter, scores in profile.iteritems():
            profile[letter] = [float(i) / self.t for i in scores]
        return profile

    # (This is a bit like a static method in other languages (but not exactly the same))
    @staticmethod
    def normalise(dist):
        # cache the sum so we don't have to calculate it repeatedly in the list
        # comprehension.
        total = sum(dist)
        return [probability / total for probability in dist]

    def get_generation_prob(self, lmer, profile):
        # the probability of the whole string is the product of the probability
        # of each letter being generated in its position over all letters
        prob = 1  # start at 1 so we can multiply
        profile_probs = [profile[lmer[i]][i] for i in xrange(self.l)]
        # they could *all* be 0!
        # so return 0 in that case.
        if all(prob == 0 for prob in profile_probs):
            return 0
        # the smallest nonzero probability
        min_prob = min(filter(lambda x: x != 0, profile_probs))
        for profile_prob in profile_probs:
            # this chance can be 0. We don't want to make the whole total 0, so
            # use a very small value to approximate 0.
            # this value is smaller than everything in the list
            if profile_prob == 0:
                prob *= min_prob * VERY_SMALL
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
        """The standard way of sampling from a distribution"""
        # take a random number between 0 and 1
        rand = random.random()
        for prob in dist:
            # and subtract a probability from it
            rand -= prob
            # if this particular probability caused the random number to drop
            # to 0 or below
            if rand <= 0:
                # then we sampled according to the distribution. We return the
                # index since we're going to use it later
                return dist.index(prob)
            # otherwise just keep going
        # if we got to here, the distribution probably wasn't notmalised
        raise MathsIsHardException, "Choosing from a distribution is hard. Did you make sure the distribution was normalised?"

    def get_motif(self, profile):
        # store the motif as a list and join it at the end since strings are
        # immutable
        motif = []
        for i in xrange(self.l):
            candidates = [(letter, profile[letter][i])
                          for letter in self.letters]
            motif.append(self.max_score(candidates)[0])
        return "".join(motif)

    # calculate the score of a profile, as defined in lectures
    def get_score(self, profile):
        # maybe these list comprehensions are getting a bit out of hand....
        # the score for a profile is the sum of the most frequent nucleotide in
        # every position in the profile's indexes
        return sum(max(index_frequencies) for index_frequencies in ((freq[i] for freq in profile.itervalues()) for i in xrange(self.l)))

    def show_motif(self, motif):
        starting_positions = self.motif_to_starting_positions[motif]
        """Show where the motif is in each sequence, in a pretty way
        eg:
            aaaaaCATaa
            ccccCATccc
            gggCATgggg
        """
        for i in xrange(self.t):
            sequence = self.DNA[i]
            starting_position = starting_positions[i]
            #lowercase if it doesn't match"".join(letter if motif[i] [sequence[starting_position: starting_position + self.l]]) + 
            print sequence[:starting_position].lower() + \
                sequence[starting_position: starting_position + self.l] + \
                sequence[starting_position + self.l:].lower()

if __name__ == "__main__":

    if len(sys.argv) < 5:
        print """
        Usage: (python|pypy) sampler.py <file> <l> <k> <n>

        file = <input file containing DNA matrix> 
        l = <length of motif to search for> 
        k = <number of identical results needed to acquire "convergence">
        n = <number of times to re-run the algorithm from a random starting point>
        """
        sys.exit(1)

    # we take the length of the motif to find as a command line argument
    l = int(sys.argv[2])
    # the number of iterations for which nothing changing will constitute
    # "convergence"
    k = int(sys.argv[3])
    n = int(sys.argv[4])
    silly = False
    if len(sys.argv) == 6:
        if sys.argv[5] == "--silly":
            silly = True
    filename = sys.argv[1]
    sampler = Sampler(filename, l, k)
    sampler.sample_random_starting_points(nsamples=n, silly=silly)

import random
import sys
class DataGenerator():
    """A class to randomly generate data for the sampler"""
    # the letters to randomly include in a sequence. We'll just assume ACTG for simplicity
    LETTERS = ['A', 'T', 'C', 'G']
    def __init__(self, t, n, motif, mutations=0):

        #the number of DNA sequences
        self.t = t
        #the length of each sequence
        self.n = n
        #the length of the motif (duh)
        self.l = len(motif)
        #the motif to include in each DNA sequence
        self.motif = motif
        #how many mutations to make in the motif
        self.mutations = mutations

    def generate_DNA_sequence(self):
        sequence = "".join(random.choice(self.LETTERS) for i in xrange(self.n - self.l))
        starting_point = random.randint(0, self.n - self.l)
        #insert the motif at a random spot
        sequence = sequence[:starting_point] + self.mutate(self.motif) + sequence[starting_point:]
        return sequence

    def mutate(self, motif):
        #keep track of where we mutated, since we can't mutate the same spot twice. That only counts as one mutation! Or even 0 if we mutate back to the orignal nucelotide!
        mutation_indicies = []
        mutated_motif = motif
        for i in xrange(self.mutations):
            #start at -1 so this definitely isn't an index in the motif
            index = -1
            #keep generating until we get something new
            while index not in mutation_indicies:
                index = random.randint(0, self.l)

            current_letter = self.mutated_motif[index]
            #choose something that isn't already at that position
            #if we wanted to be REALLY efficient we could cache this list for each letter but come on this is Python
            mutated_letter = random.choice(filter(lambda letter : letter != current_letter, self.LETTERS))
            #remember that we mutated this spot
            mutation_indicies.append(index)
            #finally mutate the letter
            mutated_motif = mutated_motif[:index] + mutated_letter + mutated_motif[index:]

        return mutated_motif

    def write_data(self, filename):
        DNA = [] 

        #if the file doesn't exist, create it, if it does, wipe it!
        open(filename, 'w').close()
        with open(filename, 'w') as file:
            for i in xrange(self.t):
                DNA.append(self.generate_DNA_sequence())

            file.write("\n".join(DNA))


if __name__ == "__main__":

    if len(sys.argv) != 6:
        print """
            Usage:
                pypy DataGenerator.py <t> <n> <motif> <number of mutations> <file path to write to>
        """
        sys.exit(1)

    t = int(sys.argv[1])
    n = int(sys.argv[2])
    motif = sys.argv[3]
    mutations = int(sys.argv[4])
    filename = sys.argv[5]

    generator = DataGenerator(t,n,motif,mutations)
    generator.write_data(filename)

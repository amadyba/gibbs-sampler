import Sampler 
import DataGenerator



def make_sample():
    t = 1000
    n = 1000
    motif = "GATTACAGATTACA"
    mutations = 0 
    generator = DataGenerator.DataGenerator(t, n, motif, mutations)
    filename = "test/%s-%s-%s-%s.txt" % (t, n, motif, mutations)
    generator.write_data(filename)

def sample_file(filename, l, k, nsamples):
    sampler = Sampler.Sampler(filename, l, k)
    motif, score = sampler.sample_random_starting_points(nsamples)
    return motif, score


def test_samples(n):
    for i in xrange(1, n, 10):
        motif, score = sample_file("test/1000-1000-GATTACAGATTACA-0.txt", 14, 10, i)
        print i, motif, score
make_sample()
test_samples(100)

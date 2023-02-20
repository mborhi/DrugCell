import numpy as np

def read_data(filename):
    samples = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            tokens = line.split("\t")
            samples.append(tokens[0])
    
    return samples

def write_mock_data(filename, samples):
    with open(filename, "w") as f:
        for sample in samples:
            type = np.random.choice([0, 1])
            line = sample + "\t" + str(type) + "\n"
            f.write(line)

if __name__ == "__main__":
    train_samples = read_data("./../data/drugcell_train.txt")
    write_mock_data("./../custom_data/mock_train_data.txt", train_samples)
    test_samples = read_data("./../data/drugcell_test.txt")
    write_mock_data("./../custom_data/mock_test_data.txt", test_samples)
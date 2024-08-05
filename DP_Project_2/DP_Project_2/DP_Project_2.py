import numpy as np

class MarkovChainAnalyzer:
    def __init__(self, transition_matrix, initial_vector):
        self.transition_matrix = np.array(transition_matrix)
        self.initial_vector = np.array(initial_vector)
        self.num_states = self.transition_matrix.shape[0]

    def classify_chain_type(self):
        if np.any(np.diag(self.transition_matrix) == 1):
            if np.all(np.diag(self.transition_matrix) == 1):
                return "Absorbing"
            else:
                return "Irreducible, Aperiodic"

        if self._check_periodicity():
            return "Aperiodic"

        if self._check_ergodicity():
            return "Ergodic"

        return "Unknown"

    def _check_periodicity(self):
        period = 1
        current_power = np.copy(self.transition_matrix)
        while not np.allclose(current_power, np.eye(self.num_states)):
            current_power = np.dot(current_power, self.transition_matrix)
            period += 1
            if period > self.num_states:
                return False
        return period % 2 == 0

    def _check_ergodicity(self):
        return np.all(np.linalg.matrix_power(self.transition_matrix, self.num_states) > 0)

    def classify_state_type(self, state_index):
        if np.all(self.transition_matrix[state_index, :] == 0) and self.transition_matrix[state_index, state_index] == 1:
            return "Absorbing"
        elif self._check_recurrent_state(state_index):
            return "Recurrent"
        else:
            return "Transient"

    def _check_recurrent_state(self, state_index):
        recurrent_states = np.where(np.linalg.matrix_power(self.transition_matrix, self.num_states)[:, state_index] > 0)[0]
        return state_index in recurrent_states

    def calculate_steady_state_probability(self):
        eigenvalues, eigenvectors = np.linalg.eig(self.transition_matrix.T)
        eigenvector = np.real(eigenvectors[:, np.isclose(eigenvalues, 1)])
        steady_state = eigenvector[:, 0] / np.sum(eigenvector[:, 0])
        return steady_state

    def calculate_state_probabilities_after_n_steps(self, n):
        state_probs = np.linalg.matrix_power(self.transition_matrix.T, n).dot(self.initial_vector)
        return state_probs

n_states = int(input("Enter the number of states: "))
transition_matrix = []
print("Enter the transition probability matrix:")
for i in range(n_states):
    row = list(map(float, input(f"Enter the transition probabilities for state {i + 1} with spaces between: ").split()))
    transition_matrix.append(row)

initial_vector = list(map(float, input("Enter the initial probability vector: ").split()))
n_steps = int(input("Enter the number of steps: "))
State_IDs=list(map(str, input(f"Enter the state ID all states with spaces between: ").split()))
markov_chain_analyzer = MarkovChainAnalyzer(transition_matrix, initial_vector)
chain_type = markov_chain_analyzer.classify_chain_type()
state_type=[]
for i in range(n_states):
    state_type.append(markov_chain_analyzer.classify_state_type(i))
steady_state_prob = markov_chain_analyzer.calculate_steady_state_probability()
n_steps_prob = markov_chain_analyzer.calculate_state_probabilities_after_n_steps(n_steps)

print("Markov Chain Type:", chain_type)
for i in range(n_states):
    print("State Type of '",State_IDs[i],"':" , state_type[i])

print("Steady State Probability:", steady_state_prob)
print(f"State Probabilities after {n_steps} steps:", n_steps_prob)

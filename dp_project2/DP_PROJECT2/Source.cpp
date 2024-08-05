#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

class MarkovChain {
public:
    MarkovChain(const std::vector<std::vector<double>>& transition_matrix, const std::vector<double>& initial_vector) {
        this->transition_matrix = Eigen::Map<const Eigen::MatrixXd>(transition_matrix[0].data(), transition_matrix.size(), transition_matrix[0].size());
        this->initial_vector = Eigen::Map<const Eigen::VectorXd>(initial_vector.data(), initial_vector.size());
        this->num_states = this->transition_matrix.rows();
    }

    std::string classify_chain() {
        if ((this->transition_matrix.diagonal().array() == 1).all()) {
            return "Absorbing";
        }
        else if ((this->transition_matrix.diagonal().array() == 1).any()) {
            return "Irreducible, Aperiodic";
        }
        else if (this->check_periodicity()) {
            return "Aperiodic";
        }
        else if (this->check_ergodicity()) {
            return "Ergodic";
        }
        else {
            return "Unknown";
        }
    }

    std::string classify_state(int state_index) {
        if ((this->transition_matrix.row(state_index).array() == 0).all() && this->transition_matrix(state_index, state_index) == 1) {
            return "Absorbing";
        }
        else if (this->check_recurrent_state(state_index)) {
            return "Recurrent";
        }
        else {
            return "Transient";
        }
    }

    Eigen::VectorXd steady_state_probability() {
        Eigen::ComplexEigenSolver<Eigen::MatrixXd> es(this->transition_matrix.transpose());
        Eigen::VectorXcd eigenvectors = es.eigenvectors().col(0);
        Eigen::VectorXd steady_state = eigenvectors.real() / eigenvectors.real().sum();
        return steady_state;
    }

    Eigen::VectorXd state_probabilities_after_n_steps(int n) {
        return (this->transition_matrix.transpose().pow(n) * this->initial_vector).real();
    }

private:
    bool check_periodicity() {
        int period = 1;
        Eigen::MatrixXd current_power = this->transition_matrix;
        while (!current_power.isApprox(Eigen::MatrixXd::Identity(this->num_states, this->num_states))) {
            current_power = current_power * this->transition_matrix;
            period++;
            if (period > this->num_states) {
                return false;
            }
        }
        return period % 2 == 0;
    }

    bool check_ergodicity() {
        return (this->transition_matrix.pow(this->num_states).array() > 0).all();
    }

    bool check_recurrent_state(int state_index) {
        Eigen::VectorXd recurrent_states = (this->transition_matrix.pow(this->num_states).col(state_index).array() > 0).cast<double>();
        return recurrent_states(state_index) > 0;
    }

    Eigen::MatrixXd transition_matrix;
    Eigen::VectorXd initial_vector;
    int num_states;
};

int main() {
    int n_states;
    std::cout << "Enter the number of states: ";
    std::cin >> n_states;

    std::vector<std::vector<double>> transition_matrix(n_states, std::vector<double>(n_states));
    std::cout << "Enter the transition probability matrix:" << std::endl;
    for (int i = 0; i < n_states; i++) {
        std::cout << "Enter the transition probabilities for state " << i + 1 << " with spaces between: ";
        for (int j = 0; j < n_states; j++) {
            std::cin >> transition_matrix[i][j];
        }
    }

    std::vector<double> initial_vector(n_states);
    std::cout << "Enter the initial probability vector: ";
    for (int i = 0; i < n_states; i++) {
        std::cin >> initial_vector[i];
    }

    int n_steps;
    std::cout << "Enter the number of steps: ";
    std::cin >> n_steps;

    MarkovChain chain(transition_matrix, initial_vector);
    std::string chain_type = chain.classify_chain();
    std::string state_type = chain.classify_state(0);
    Eigen::VectorXd steady_state = chain.steady_state_probability();
    Eigen::VectorXd state_probs_n_steps = chain.state_probabilities_after_n_steps(n_steps);

    std::cout << "Markov Chain Type: " << chain_type << std::endl;
    std::cout << "State Type: " << state_type << std::endl;
    std::cout << "Steady State Probability: " << steady_state.transpose() << std::endl;
    std::cout << "State Probabilities after " << n_steps << " steps: " << state_probs_n_steps.transpose() << std::endl;

    return 0;
}


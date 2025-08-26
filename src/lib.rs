use std::vec;

use complex_rs::Complex;
use rain_linalg::Vector;
use num_traits::Float;
use rain_linalg::Matrix;
use rand::{distr::Distribution, rng};

#[derive(Debug)]
pub struct QuantumRegister<T: Float> {
    num_qubits: usize,
    state_vector: Vector<T>,
}

impl<T: Float> QuantumRegister<T> {
    pub fn new(num_qubits: usize) -> Self {
        let dim = 1 << num_qubits;
        let mut elements = vec![Complex::zero(); dim];
        elements[0] = Complex::one();
        Self {
            num_qubits,
            state_vector: Vector::new(elements),
        }
    }

    pub fn state_vector(&self) -> &Vector<T> {
        &self.state_vector
    }

    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }



    pub fn apply_gate(&mut self, targets: &[usize], gate: &QuantumGate<T>) {
    let num_targets = targets.len();
    let gate_dim = 1 << num_targets;
    assert_eq!(gate.matrix.rows(), gate_dim, "Gate dimension does not match number of target qubits.");

    let num_other_qubits = self.num_qubits - num_targets;
    let mut other_qubits = Vec::with_capacity(num_other_qubits);
    for i in 0..self.num_qubits {
        if !targets.contains(&i) {
            other_qubits.push(i);
        }
    }

    let mut new_state_elements = self.state_vector.elements().to_vec();

    for i in 0..(1 << num_other_qubits) {
        let mut sub_vector_elements = Vec::with_capacity(gate_dim);
        let mut original_indices = Vec::with_capacity(gate_dim);

        for j in 0..gate_dim {
    
            let mut index = 0;
            let mut other_pos = 0;
            for k in 0..self.num_qubits {
                if let Some(pos) = targets.iter().position(|&t| t == k) {
                    let bit_pos = num_targets - 1 - pos;
                    if (j >> bit_pos) & 1 == 1 {
                        index |= 1 << k;
                    }
                } else {
                    if (i >> other_pos) & 1 == 1 {
                        index |= 1 << k;
                    }
                    other_pos += 1;
                }
            }
            
            sub_vector_elements.push(self.state_vector.get(index));
            original_indices.push(index);
        }

        let sub_vector = Vector::new(sub_vector_elements);
        let transformed_sub_vector = gate.matrix.clone() * sub_vector;

        for j in 0..gate_dim {
            new_state_elements[original_indices[j]] = transformed_sub_vector.get(j);
        }
    }
    
    self.state_vector = Vector::new(new_state_elements);
}

    pub fn measure(&self) -> usize {
        let probabilities: Vec<f64> = self.state_vector.elements().iter().map(|c| c.magnitude_squared().to_f64().unwrap()).collect();

        let mut rng = rng();
        let dist = 
        rand::distr::weighted::WeightedIndex::new(&probabilities).unwrap();
        dist.sample(&mut rng)

    }
}

#[derive(Debug, Clone)]
pub struct QuantumGate<T: Float> {
    pub name: String,
    matrix: Matrix<T>,
}

// matrix must be unitary matrix!
impl<T: Float> QuantumGate<T> {
    pub fn new(name:String,matrix: Matrix<T>) -> Self {
    
        Self { name, matrix }
    }


    /// Identity Gate
    pub fn i() -> Self {
        Self::new("I".to_string(), Matrix::new(2,2, vec![
            Complex::one(), Complex::zero(),
            Complex::zero(), Complex::one(),
        ]))
    }

    /// NOT Gate ( Pauli - X )
    pub fn x() -> Self {
        Self::new("X".to_string(), Matrix::new(2,2, vec![
            Complex::zero(), Complex::one(),
            Complex::one(),  Complex::zero(),
        ]))
    }

    /// Pauli - Y Gate
    pub fn y() -> Self {
        Self::new("Y".to_string(), Matrix::new(2,2, vec![
            Complex::zero(), -Complex::i(),
            Complex::i(),    Complex::zero(),
        ]))
    }

    /// Pauli - Z Gate
    pub fn z() -> Self {
        Self::new("Z".to_string(), Matrix::new(2,2, vec![
            Complex::one(),     Complex::zero(),
            Complex::zero(),    -Complex::one(),
        ]))
    }

    /// Hadmard Gate
    pub fn h() -> Self {
        let factor = T::one() / T::from(2.0).unwrap().sqrt();
        Self::new("H".to_string(), Matrix::new(2,2, vec![
            Complex::one() * factor, Complex::one() * factor,
            Complex::one() * factor, -Complex::one() * factor,
        ]))
    }

    /// Phase Gate
    pub fn s() -> Self {
        Self::new("S".to_string(), Matrix::new(2,2, vec![
            Complex::one(),     Complex::zero(),
            Complex::zero(),    Complex::i(),
        ]))
    }

    /// T gate
    pub fn t() -> Self {
        let pi = T::from(std::f64::consts::PI).unwrap();
        let angle = pi / T::from(4.0).unwrap();
        Self::new("T".to_string(), Matrix::new(2,2, vec![
            Complex::one(),     Complex::zero(),
            Complex::zero(),    Complex::new(angle.cos(), angle.sin()),
        ]))
    }

    /// Controlled-NOT Gate
    pub fn cnot() -> Self {
        Self::new("CNOT".to_string(), Matrix::new(4,4, vec![
            Complex::one(),     Complex::zero(),    Complex::zero(),    Complex::zero(),
            Complex::zero(),    Complex::one(),     Complex::zero(),    Complex::zero(),
            Complex::zero(),    Complex::zero(),    Complex::zero(),    Complex::one(), 
            Complex::zero(),    Complex::zero(),    Complex::one(),     Complex::zero(),
        ]))
    }
}

pub struct QuantumCircuit<T: Float> {
    steps: Vec<(QuantumGate<T>, Vec<usize>)>,
}

impl<T: Float> QuantumCircuit<T> {
    pub fn new() -> Self {
        Self { steps: Vec::new()}
    }

    pub fn add_gate(&mut self, targets: &[usize], gate:QuantumGate<T>) {
        self.steps.push((gate, targets.to_vec()));
    }

    pub fn run(&self, register: &mut QuantumRegister<T>) {
        for (gate, targets) in &self.steps {
            register.apply_gate(targets, gate);
        }
    }
}




#[cfg(test)]
mod tests
{
    use super::*;
    use rain_linalg::Vector;
    use super::QuantumCircuit;
    #[test]
    fn test_new_quantum_register() {
        let reg1 = QuantumRegister::<f64>::new(1);
        let expected1 = Vector::new(vec![Complex::one(), Complex::zero()]);
        assert_eq!(reg1.state_vector(), &expected1);
        assert_eq!(reg1.num_qubits(), 1);

        let reg2 = QuantumRegister::<f64>::new(2);
        let expected2 = Vector::new(vec![
            Complex::one(),
            Complex::zero(),
            Complex::zero(),
            Complex::zero(),
        ]);
        assert_eq!(reg2.state_vector(), &expected2);
        assert_eq!(reg2.num_qubits(), 2);
    }
        const TOLERANCE: f64 = 1e-10;

    fn assert_complex_eq(a: Complex<f64>, b: Complex<f64>) {
        assert!((a.re - b.re).abs() < TOLERANCE);
        assert!((a.im - b.im).abs() < TOLERANCE);
    }
    
    fn assert_vector_eq(a: &Vector<f64>, b: &Vector<f64>) {
        assert_eq!(a.dim(), b.dim());
        for (ca, cb) in a.elements().iter().zip(b.elements().iter()) {
            assert_complex_eq(*ca, *cb);
        }
    }

    // ... 이전 QuantumRegister 테스트 ...

    #[test]
    fn test_gate_creation() {
        let gate_x = QuantumGate::<f64>::x();
        assert_eq!(gate_x.name, "X");
        assert_eq!(gate_x.matrix.get(0, 1), Complex::one());
    }

    #[test]
    fn test_hadamard_on_zero() {
        let h_gate = QuantumGate::<f64>::h();
        let state_zero = Vector::new(vec![Complex::one(), Complex::zero()]);
        
        // H|0> = |+> = 1/sqrt(2) * [1, 1]
        let result_vec = h_gate.matrix * state_zero;
        
        let factor = 1.0 / 2.0_f64.sqrt();
        let expected_vec = Vector::new(vec![
            Complex::new(factor, 0.0),
            Complex::new(factor, 0.0)
        ]);

        assert_vector_eq(&result_vec, &expected_vec);
    }

    #[test]
    fn test_cnot_on_10() {
        let cnot_gate = QuantumGate::<f64>::cnot();
        // |10> 상태 벡터 = [0, 0, 1, 0]
        let state_10 = Vector::new(vec![
            Complex::zero(), Complex::zero(), Complex::one(), Complex::zero()
        ]);

        // CNOT|10> = |11>
        let result_vec = cnot_gate.matrix * state_10;
        
        // |11> 상태 벡터 = [0, 0, 0, 1]
        let expected_vec = Vector::new(vec![
            Complex::zero(), Complex::zero(), Complex::zero(), Complex::one()
        ]);

        assert_vector_eq(&result_vec, &expected_vec);
    }
    

        #[test]
    fn test_bell_state_creation_fully() {
        let mut reg = QuantumRegister::<f64>::new(2);
        reg.apply_gate(&[0], &QuantumGate::h());

        // CNOT (control=0, target=1) -> targets를 [0, 1] 순서로 전달
        reg.apply_gate(&[0, 1], &QuantumGate::cnot());

        let factor = 1.0 / 2.0_f64.sqrt();
        let expected_vec = Vector::new(vec![
            Complex::new(factor, 0.0),
            Complex::zero(),
            Complex::zero(),
            Complex::new(factor, 0.0),
        ]);
        assert_vector_eq(reg.state_vector(), &expected_vec);
    }

    #[test]
    fn test_ghz_state_creation() {
        let mut reg = QuantumRegister::<f64>::new(3);
        reg.apply_gate(&[0], &QuantumGate::h());

        // CNOT (control=0, target=1) -> targets: [0, 1]
        reg.apply_gate(&[0, 1], &QuantumGate::cnot());

        // CNOT (control=0, target=2) -> targets: [0, 2]
        reg.apply_gate(&[0, 2], &QuantumGate::cnot());

        let factor = 1.0 / 2.0_f64.sqrt();
        let mut expected_elements = vec![Complex::zero(); 8];
        expected_elements[0] = Complex::new(factor, 0.0);
        expected_elements[7] = Complex::new(factor, 0.0);
        let expected_vec = Vector::new(expected_elements);
        assert_vector_eq(reg.state_vector(), &expected_vec);

    }

    #[test]
fn test_quantum_circuit_run() {
    
    let mut bell_circuit = QuantumCircuit::<f64>::new();
    
    bell_circuit.add_gate(&[0], QuantumGate::h());
    
    
    bell_circuit.add_gate(&[0, 1], QuantumGate::cnot());

    let mut register = QuantumRegister::<f64>::new(2);

    
    bell_circuit.run(&mut register);

    
    let factor = 1.0 / 2.0_f64.sqrt();
    let expected_vec = Vector::new(vec![
        Complex::new(factor, 0.0),
        Complex::zero(),
        Complex::zero(),
        Complex::new(factor, 0.0),
    ]);

    assert_vector_eq(register.state_vector(), &expected_vec);
}
}


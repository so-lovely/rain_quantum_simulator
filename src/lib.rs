use std::vec;

use complex_rs::Complex;
use rain_linalg::Vector;
use num_traits::Float;
use rain_linalg::Matrix;

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
}

#[derive(debug)]
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


#[cfg(test)]
mod tests
{
    use super::*;
    use rain_linalg::Vector;
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
        for (ca, cb) in a.elements.iter().zip(b.elements.iter()) {
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
}

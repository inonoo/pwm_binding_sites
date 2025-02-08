import math
import numpy as np

class PWM():
  def __init__(self, path: str):
    # Initialize all parameters
    self.mat = self._read_pwm(path)

  def find(self, seq: str, threshold: float = 1.0e-05) -> list[tuple[int, float]]:
    """From a sequence, find binding sites and their scores."""
    out = []
    pwm_length = len(self.mat[list(self.mat.keys())[0]])
    for i in range(len(seq) - pwm_length):  # don't run off the ends
      sub_seq = seq[i:i + pwm_length]
      sub_score = self.score(sub_seq)
      if sub_score > threshold:
        out.append((i, sub_score))
    return out

  def score(self, seq: str) -> float:
    """Return the matching score/probability for a single position. In this case
    the sequence length must be IDENTICAL to the PWM length."""
    assert len(seq) == len(self.mat[list(self.mat.keys())[0]])
    seq = seq.lower()

    score = 0
    for i, letter in enumerate(seq):
      score += self.mat[letter][i]

    score_rev = 0
    for i, letter in enumerate(self._invert_seq(seq)):
      score_rev += self.mat[letter][i]

    score = score if score > score_rev else score_rev
    return math.exp(score)

  def _read_pwm(self, path: str):
    """Read in a PWM file. Maybe a numpy matrix would be better, but right now
    we're going to use a dict."""
    # This is inefficient. If you want to efficient, consider turning a's into
    # 0s, c's into 1s, and then do everything in numpy matrix by position.
    out = {}
    with open(path, 'r') as file:
      for i, line in enumerate(file.readlines()):
        if i == 0:
          alphabet = line[10:].lower()
          for letter in alphabet:
            out[letter] = []
        elif i > 1:
          line = [int(v) for v in line.strip().split(' ') if len(v) > 0]
          for letter, prob in zip(alphabet, line):
            out[letter].append(prob)
    return out

  def read_as_numpy(path: str):
      """
      Reads a file and returns a NumPy matrix.

      Args:
          path (str): The file path.

      Returns:
          np.ndarray: A NumPy matrix.
      """
      with open(path, 'r') as file:
          lines = file.readlines()

      # Extract the alphabet from the first line
      alphabet = lines[0].split('=')[1].strip()

      # The numeric matrix
      # (skip the first two lines)
      matrix = np.array([list(map(int, line.split())) for line in lines[2:]])

      return matrix

  def _invert_seq(self, seq: str) -> str:
    """Convert As to Ts and Cs to Gs and reverse."""
    # If this was in numpy, we could just do 3 - the number and reverse order
    out = ''
    for letter in seq[::-1]:  # reverse a list
      if letter == 'a':
        out += 't'
      elif letter == 'c':
        out += 'g'
      elif letter == 'g':
        out += 'c'
      elif letter == 't':
        out += 'a'
      else:
        raise ValueError('I cannot handle these crazy letters')
    return out

with open('ebv.txt', 'r') as file:
  seq = file.read().strip()

pwm = PWM('pwm.txt')
print(pwm.find(seq))
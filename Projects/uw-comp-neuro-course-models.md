# University of Washington Computational Neuroscience Course Models

## Short Description

This project collects selected Python models I built while working through the University of Washington online computational neuroscience course.

Instead of showing every small test file, I grouped the strongest work into four themes:

- neural coding and population decoding
- membrane dynamics and spiking models
- synaptic input and spike-rate adaptation
- eigenvectors and learning from input structure

I cannot share the datasets used in these models so they will not run. They are here to demonstrate coding and problem-solving abilites and computational neuroscience knowledge.

---

## 1. Neural Coding and Population Decoding

This section combines two related ideas:

- checking whether neural firing variability looks Poisson-like
- using a population vector to decode a stimulus direction

Together, these scripts show how neural responses can be analyzed statistically and then used to estimate information about a stimulus.

### Poisson Variability Check

```python
import pickle
import matplotlib.pyplot as plt
import numpy as np

with open('tuning_3.4.pickle', 'rb') as f:
    data = pickle.load(f)

all_means = []
all_vars = []

for neuron in ['neuron1', 'neuron2', 'neuron3', 'neuron4']:
    mean_fire = np.mean(data[neuron], axis=0)
    var_fire = np.var(data[neuron], axis=0)
    all_means.extend(mean_fire)
    all_vars.extend(var_fire)

    plt.scatter(mean_fire, var_fire, label=neuron)

max_axis = max(max(all_means), max(all_vars))
plt.plot([0, max_axis], [0, max_axis], '--')
plt.xlabel('Mean firing rate')
plt.ylabel('Variance')
plt.title('Poisson check: variance vs mean')
plt.legend()
plt.show()
```
![](Figures/poisson_check.png)

### Population Vector Decoding

```python
import pickle
import numpy as np

with open('pop_coding_3.4.pickle', 'rb') as f:
    datap = pickle.load(f)

with open('tuning_3.4.pickle', 'rb') as f:
    data = pickle.load(f)

rmax1 = np.max(np.mean(data['neuron1'], axis=0))
rmax2 = np.max(np.mean(data['neuron2'], axis=0))
rmax3 = np.max(np.mean(data['neuron3'], axis=0))
rmax4 = np.max(np.mean(data['neuron4'], axis=0))

max_rates = np.array([rmax1, rmax2, rmax3, rmax4])
if not np.all(np.isfinite(max_rates) & (max_rates > 0)):
    raise ValueError("Maximum firing rates must be finite and positive.")

r1 = np.mean(datap['r1'])
r2 = np.mean(datap['r2'])
r3 = np.mean(datap['r3'])
r4 = np.mean(datap['r4'])

v = (r1 / rmax1) * datap['c1'] \
    + (r2 / rmax2) * datap['c2'] \
    + (r3 / rmax3) * datap['c3'] \
    + (r4 / rmax4) * datap['c4']

if not np.all(np.isfinite(v)) or np.all(v == 0):
    raise ValueError("Cannot decode a direction from a zero or non-finite vector.")

angle = np.degrees(np.arctan2(v[1], v[0]))
angle = angle % 360

print("population vector:", v)
print("angle:", round(angle))
```

---

## 2. Membrane Dynamics and Spiking Models

This section combines passive membrane charging, integrate-and-fire dynamics, and noisy spike timing.

The core idea is that a neuron's voltage changes over time according to input current, leak, threshold, reset, and noise.

### Passive Membrane Time Constant

```python
import numpy as np
import matplotlib.pyplot as plt

I = 10  # nA
C = 0.1  # nF
R = 100  # M ohms

tau_theoretical = R*C

print("C =", C, "nF")
print("R =", R, "M ohms")
print("tau =", tau_theoretical, "ms")
print('(Theoretical)')

tstop = 150
V_inf = I*R
tau_numerical = None
h = 0.2
V = 0
V_trace = [V]
time_points = [0]

for t in np.arange(h, tstop, h):
    V = V + h*(- (V/(R*C)) + (I/C))

    
    if tau_numerical is None and V > 0.6321*V_inf:
        tau_numerical = t
        print("tau =", round(tau_numerical, 3), "ms")
        print('(Numerical estimate)')

    if t >= 0.6*tstop:
        I = 0

    V_trace += [V]
    time_points += [t]

plt.plot(time_points, V_trace)
plt.xlabel("Time (ms)")
plt.ylabel("Voltage relative to rest (mV)")
plt.title("Passive membrane charging")
plt.show()
```

### Integrate-and-Fire Neuron

```python
import numpy as np
import matplotlib.pyplot as plt

I = 10
C = 1
R = 40

V = 0
tstop = 200
abs_ref = 5
ref = 0
V_trace = [V]
time_points = [0]
V_th = 10

for t in range(tstop):

    if not ref:
        V = V - (V/(R*C)) + (I/C)
    else:
        ref -= 1
        V = 0.2 * V_th

    if V > V_th:
        V = 50
        ref = abs_ref

    V_trace += [V]
    time_points += [t + 1]

plt.plot(time_points, V_trace)
plt.xlabel("Time (ms)")
plt.ylabel("Voltage")
plt.title("Integrate-and-fire neuron")
plt.show()
```

### Noise and Interspike Intervals

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0)

C = 1
R = 40
tstop = 10000
abs_ref = 5
V_th = 10
baseline_I = 1

noise_values = np.arange(0, 5.5, 0.5)

for noiseamp in noise_values:

    V = 0
    ref = 0
    spiketimes = []
    I = baseline_I + noiseamp * np.random.normal(0, 1, tstop)

    for t in range(tstop):

        if ref == 0:
            V = V - V / (R * C) + I[t] / C
        else:
            ref -= 1
            V = 0.2 * V_th

        if V > V_th:
            V = 50
            ref = abs_ref
            spiketimes.append(t + 1)

    isi = np.diff(spiketimes)

    if len(isi) >= 2:
        if np.ptp(isi) == 0:
            plt.axvline(isi[0], linestyle="--", color="black",
                        label=f"Noise = {noiseamp:.1f}: {isi[0]} ms (deterministic)")
            continue
        bandwidth = 2
        x = np.linspace(max(0, isi.min() - 4 * bandwidth),
                        isi.max() + 4 * bandwidth, 600)
        y = np.zeros_like(x)

        for value in isi:
            y += np.exp(-0.5 * ((x - value) / bandwidth) ** 2)

        y = y / (len(isi) * bandwidth * np.sqrt(2 * np.pi))
        plt.plot(x, y, label=f"Noise = {noiseamp:.1f}")

plt.xlabel("Interspike interval (ms)")
plt.ylabel("Probability density")
plt.title("ISI distributions across noise amplitudes")
plt.legend(fontsize=7)
plt.tight_layout()
from pathlib import Path
Path("Figures").mkdir(exist_ok=True)
plt.savefig("Figures/ISI_distributions.png", dpi=160)
plt.show()
```
![ISI distributions](Figures/ISI_distributions.png)

---

## 3. Synaptic Input and Spike-Rate Adaptation

This model is one of the strongest pieces in the folder because it goes beyond a single current input.

It uses:

- a random input spike train
- alpha-function synaptic conductance
- excitatory synaptic current
- refractory periods
- spike-rate adaptation
- spike count as an output measure

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0)
h = 1.
t_max = 200
tstop = int(t_max/h)

input_rate_hz = 100
thr = 1 - input_rate_hz * h / 1000  # Bernoulli probability per step
spike_train = np.random.rand(tstop) > thr

t_a = 100
t_peak_values = np.arange(0.5, 10.5, 0.5)
g_peak = 0.05
t_vec = np.arange(0, t_a + h, h)

C = 0.5
R = 40
G_inc = 1.0  # fixed conductance increment per output spike
tau_ad = 2

E_leak = -60
E_syn = 0
V_th = -40
V_spike = 50
ref_max = int(4/h)

fig, axs = plt.subplots(2, 1)
axs[0].plot(np.arange(0, t_max, h), spike_train)
axs[0].set_title('Input spike train')

spike_counts = []

for t_peak in t_peak_values:

    const = g_peak / (t_peak * np.exp(-1))
    alpha_func = const * t_vec * np.exp(-t_vec / t_peak)

    V = E_leak
    g_ad = 0
    t_list = np.array([], dtype=int)
    spike_count = 0
    ref = 0

    V_trace = [V]
    t_trace = [0]

    for t in range(tstop):

        if spike_train[t]:
            t_list = np.concatenate([t_list, [0]])

        g_syn = np.sum(alpha_func[t_list])
        I_syn = g_syn*(E_syn - V)

        if t_list.size > 0:
            t_list = t_list + 1
            t_list = t_list[t_list < len(alpha_func)]

        if not ref:
            V = V + h*(-((V-E_leak)*(1+R*g_ad)/(R*C)) + (I_syn/C))
            g_ad = g_ad + h*(-g_ad/tau_ad)
        else:
            ref -= 1
            V = V_th - 10
            g_ad = g_ad + h*(-g_ad/tau_ad)

        if (V > V_th) and not ref:
            V = V_spike
            spike_count += 1
            ref = ref_max
            g_ad = g_ad + G_inc

        V_trace += [V]
        t_trace += [(t + 1)*h]

    spike_counts.append(spike_count)

plt.figure()
plt.plot(t_peak_values, spike_counts, marker="o")
plt.xlabel("Synaptic time to peak")
plt.ylabel("Spike count")
plt.title("Spike count across synaptic time constants")

axs[1].plot(t_trace, V_trace)
axs[1].set_title('Output voltage trace')
plt.show()

print("Spike Count:", spike_counts)
```

![Spike Count Vs. Synaptic time to peak](Figures/Spike_count.png)

---

## 4. Eigenvectors and Learning From Input Structure

This section combines the input-correlation eigenvector script with the centered-data learning script.

The shared idea is that a system can learn or identify important directions in input data.

```python
import pickle
import matplotlib.pyplot as plt
import numpy as np

with open('c10p1.pickle', 'rb') as f:
    data = pickle.load(f)

points = data["c10p1"]

mean_point = np.mean(points, axis=0)
u = points - mean_point

print("Mean after centering:", np.mean(u, axis=0))

plt.scatter(u[:, 0], u[:, 1])
plt.axhline(0)
plt.axvline(0)
plt.xlabel("Centred x")
plt.ylabel("Centred y")
plt.title("Mean-centred input data")
plt.show()

eta = 1
dt = 0.01
num_iterations = 100_000

rng = np.random.default_rng(0)
w = rng.random(2)
w = w / np.linalg.norm(w)
num_points = len(u)

for t in range(num_iterations):
    index = t % num_points
    current_u = u[index]
    v = current_u @ w

    delta_w = dt*eta * (v * current_u - (v**2) * w)
    w = w + delta_w

print("Final weight vector:", w)
print("Length of final vector:", np.linalg.norm(w))

Q = np.cov(u.T)

vals, vecs = np.linalg.eig(Q)
principal_eigenvector = vecs[:, np.argmax(vals)]

print("Principal eigenvector:", principal_eigenvector)
print("Length of principal eigenvector:", np.linalg.norm(principal_eigenvector))
```
![Mean-centred input data](Figures/mean_centred.png)

---

# xas-simulator

[Nexpy](https://github.com/nexpy/nexpy) plugin to simulate XAS spectra

It uses [Quanty](https://www.quanty.org/) to do multiplet calculations under different approximations.

### Installation
Requries: **Quanty**, **python >3.10**, *nexpy*, *HdfMap*, *tabulate*

**Quanty** must be downloaded from [www.quanty.org](https://www.quanty.org/), which requires an email sign-on.

```bash
$ python -m pip install nexpy hdfmap
$ git clone https://github.com/DanPorter/xas-simulator.git
$ cd xas-simulator
$ python -m pip install .
```

### Run
```bash
$ nexpy
```
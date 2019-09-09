if (utils.scm_checkout()) return

withCredentials([string(
    credentialsId: 'jwst-codecov',
    variable: 'codecov_token')]) {

env_vars = [
    "CRDS_SERVER_URL=https://jwst-crds.stsci.edu",
    "CRDS_PATH=./crds_cache",
]

// pip related setup for local index, not used currently
def pip_index = "https://bytesalad.stsci.edu/artifactory/api/pypi/datb-pypi-virtual/simple"
def pip_install_args = "--index-url ${pip_index} --progress-bar=off"


// Generate distributions build
bc0 = new BuildConfig()
bc0.nodetype = 'linux'
bc0.name = 'wheel-sdist'
bc0.conda_ver = '4.6.14'
bc0.conda_packages = [
    "python=3.6",
]
bc0.build_cmds = [
    "pip install numpy",
    "pip wheel .",
    "python setup.py sdist",
]

// Generate pip build/test with released upstream dependencies
bc1 = utils.copy(bc0)
bc1.name = "stable-deps"
bc1.env_vars = env_vars
bc1.build_cmds = [
    "pip install -e .[test]",
]
bc1.test_cmds = [
    "pytest --cov=./ -r sx --junitxml=results.xml",
    "codecov --token=${codecov_token}"
]

// Generate conda-free build with python 3.7
bc2 = new BuildConfig()
bc2.nodetype = 'python3.7'
bc2.name = 'conda-free'
bc2.env_vars = env_vars
bc2.build_cmds = [
    "pip install -e .[test]",
]
bc2.test_cmds = [
    "pytest --cov=./ -r sx --junitxml=results.xml",
    "codecov --token=${codecov_token}"
]


// Audit code with bandit
bc3 = new BuildConfig()
bc3.nodetype = 'python3.7'
bc3.name = 'security-audit'
bc3.build_cmds = [
    "pip install -e .[test] bandit",
]
bc3.test_cmds = [
    "bandit --exclude 'jwst/*test*' -r jwst scripts",
]

utils.run([bc0, bc1, bc2, bc3])
}  // withCredentials

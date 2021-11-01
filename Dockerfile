FROM python:3

WORKDIR /home/arthur/ms2rescore-immunopeptidomics-manuscript

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY immuno_ms2rescore_tools/. ./immuno_ms2rescore_tools/
COPY setup.py ./

run pip install -e .


#include <Python.h>

extern void mixtures(double *data, double *weights, double *Mus,
                     double *Covs, double *priors, int *pn,
                     int *pp, double *pdmover2, int *pK,
                     int *pKmin, double *ptol, double *pmindl,
                     int *pcountf, double *dl, double *logliks,
                     int *kappas, int *trans1, int *trans2,
                     int *lives, double *bMus, double *bCovs,
                     double *bpriors, int *pitmax, int *pverbose);

static PyObject* Py_call_mixtures(PyObject *self, PyObject *args) {
  PyArrayObject *data, *weights, *Mus, *Covs, *priors,
    *kappas, *trans1, *trans2, *lives, *logliks, *bMus, *bCovs, *bpriors;
  int *pn, *pp, *pK, *pKmin, *pcountf, *pitmax, *pverbose;
  double *pdmover2, *ptol, *pmindl, *dl;

  if (!PyArg_ParseTuple(args,
                        "OOOOOiidiiddidOOOOOOOOii",
                        &data, &weights, &Mus, &Covs,
                        &priors, pn, pp, pdmover2, pK,
                        pKmin, ptol, pmindl, pcountf, dl,
                        &logliks, &kappas, &trans1, &trans2,
                        &lives, &bMus, &bCovs, &bpriors, pitmax,
                        pverbose
                        ))
    return NULL;

  mixtures((double *)data->data,
           (double *)weights->data,
           (double *)Mus->data,
           (double *)Covs->data,
           (double *)priors->data,
           pn,
           pp,
           pdmover2,
           pK,
           pKmin,
           ptol,
           pmindl,
           pcountf,
           dl,
           (double *)logliks->data,
           (int *)kappas->data,
           (int *)trans1->data,
           (int *)trans2->data,
           (int *)lives->data,
           (double *)bMus->data,
           (double *)bCovs->data,
           (double *)bpriors->data,
           pitmax,
           pverbose);

  return Py_BuildValue(""); //何かしら返す必要があるので。無いとsegmentation fault
}

static PyMethodDef MixturesMethods[] =
  {
   {"mixtures", Py_call_mixtures, METH_VARARGS},
   {NULL},
  };

static struct PyModuleDef mixturesmodule =
  {
   PyModuleDef_HEAD_INIT,
   "mixtures",
   NULL,
   -1,
   MixturesMethods
  };

PyMODINIT_FUNC PyInit_mixtures(void)
{
  return PyModule_Create(&mixturesmodule);
}

/* HACK: This is a workaround for enforcing platform-specific wheels in Python 3.12 */
#include <stdio.h>
#include <Python.h>

void dummy() {
    printf("This is a dummy function to force a platform-specific wheel.\n");
}

/* NOTE: inspired by https://gist.github.com/douglas-larocca/099bf7460d853abb7c17 */
PyMODINIT_FUNC PyInit__native_build(void)
{

    PyObject *module;
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_chi2",
        module_docstring,
        -1,
        module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
    module = PyModule_Create(&moduledef);
    if (!module) return NULL;

    return module;
}
devtools::load_all()
Rprof("prof.out", interval = 0.01)

# Run a representative heavy task, e.g., the tests
devtools::test()

Rprof(NULL)
summary_res <- summaryRprof("prof.out")
print(head(summary_res$by.total, 20))
print(head(summary_res$by.self, 20))

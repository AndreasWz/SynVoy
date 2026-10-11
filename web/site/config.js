// Where the search service of this site is.
//
//   apiBase: ""                      the service runs on this site (the usual case)
//   apiBase: "https://api.example"   the service runs elsewhere (list this site under
//                                    "cors_origins" in the service's settings)
//   apiBase: null                    there is no service: the pages show the examples only
//
// This is the only file that differs between those three ways of hosting the pages.
window.SYNVOY = { apiBase: "" };

process ESTIMATE_PARAMS {
    tag "llm_params"
    label 'process_low'

    input:
    path resolved_query_json
    val home_species
    val target_species
    path query_fasta

    output:
    path "estimated_params.json", emit: params_json

    script:
    def target_arg = target_species ? "--target_species \"${target_species}\"" : ""
    """
    # Build biological context from resolved query + species info
    ${projectDir}/bin/build_llm_context.py \\
        --resolved_query ${resolved_query_json} \\
        --home_species "${home_species}" \\
        ${target_arg} \\
        --fasta ${query_fasta} \\
        --output context.json

    # Estimate parameters using LLM (Ollama → Google Cloud → Heuristic)
    ${projectDir}/bin/llm_param_advisor.py \\
        --context context.json \\
        --model ${params.llm_model} \\
        --ollama_url ${params.ollama_url} \\
        --ollama_timeout ${params.ollama_timeout} \\
        --google_api_key "\${GOOGLE_API_KEY:-}" \\
        --output estimated_params.json
    """
}

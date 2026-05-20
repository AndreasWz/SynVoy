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
    def api_key_arg = params.llm_api_key ? "--api_key \"${params.llm_api_key}\"" : ""
    def base_url_arg = params.llm_api_base_url ? "--api_base_url \"${params.llm_api_base_url}\"" : ""
    def model_arg = params.llm_model ? "--model \"${params.llm_model}\"" : ""
    """
    # Build biological context from resolved query + species info
    ${projectDir}/bin/build_llm_context.py \\
        --resolved_query ${resolved_query_json} \\
        --home_species "${home_species}" \\
        ${target_arg} \\
        --fasta ${query_fasta} \\
        --output context.json

    # Estimate parameters using LLM API (falls back to heuristics if no key)
    LLM_API_KEY="\${LLM_API_KEY:-\${GOOGLE_API_KEY:-\${OPENAI_API_KEY:-}}}" \\
    ${projectDir}/bin/llm_param_advisor.py \\
        --context context.json \\
        --provider ${params.llm_provider} \\
        ${api_key_arg} \\
        ${base_url_arg} \\
        ${model_arg} \\
        --output estimated_params.json
    """
}

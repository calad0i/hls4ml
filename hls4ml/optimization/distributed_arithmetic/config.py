class _global_config:
    trace_depth = False
    fuse_associative_ops = True
    backend = 'da4ml'
    order_metrics = ()
    enabled = True
    dsp_offload_thres = -1
    minimal_latency_compile = False
    use_ternary = False
    allow_split = True
    enable_pixel_unroll = False
    inp_n_max = -1
    out_n_max = -1
    n_beams = 1


class VariableOverrideContextManager:
    target = ''

    def __init__(self, value=None):
        self.value = value

    def __enter__(self):
        self.original = getattr(_global_config, self.target)
        setattr(_global_config, self.target, self.value)

    def __exit__(self, exc_type, exc_value, traceback):
        setattr(_global_config, self.target, self.original)

    def __new__(cls, *args, **kwargs):
        if cls is VariableOverrideContextManager:
            raise TypeError('VariableOverrideContextManager should not be instantiated directly')
        return super().__new__(cls)


def compiler_config(
    enabled=True,
    dsp_offload_thres: int = -1,
    minimal_latency: bool = False,
    use_ternary: bool = False,
    allow_split: bool = True,
    flatten_conv: bool = False,
    backend: str = 'da4ml',
    inp_n_max: int = -1,
    out_n_max: int = -1,
    n_beams: int = 1,
):

    _global_config.enabled = enabled
    _global_config.dsp_offload_thres = dsp_offload_thres
    _global_config.minimal_latency_compile = minimal_latency
    _global_config.inp_n_max = inp_n_max
    _global_config.out_n_max = out_n_max
    _global_config.use_ternary = use_ternary
    _global_config.allow_split = allow_split
    _global_config.enable_pixel_unroll = flatten_conv
    _global_config.backend = backend
    _global_config.n_beams = n_beams
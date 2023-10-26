import re
from typing import Optional

from hls4ml.model.types import FixedPrecisionType
from hls4ml.backends.fpga.fpga_types import NamedType
from hls4ml.model.layers import Layer, register_layer
from hls4ml.model.optimizer import register_pass, OptimizerPass

re_purge_prefix = re.compile(r'(?<!\w)(?:ap_|ac_)', re.IGNORECASE)
re_parse_fixed = re.compile(r'\s*(u?)fixed<([\w,]+)>\s*', re.IGNORECASE)


class FixedPointQuantizer(Layer):
    def initialize(self):
        inp = self.get_input_variable()
        shape = inp.shape
        dims = inp.dim_names
        self.add_output_variable(shape, dims)
        self.set_attr('n_in', self.get_input_variable().size())
        self.overrides = self.attributes['overrides']
        self.removable = self.attributes['removable']
        self.SAT, self.RND = self.attributes['SAT'], self.attributes['RND']
        self.mask_kbi = self.attributes.get('mask_kbi', None)


def to_hls4ml_fixed(fixed: str):
    matched = re_parse_fixed.match(re_purge_prefix.sub('', fixed))
    assert matched is not None, f'Cannot parse {fixed}'
    signed = matched.group(1) != 'u'
    b, i, *args = matched.group(2).split(',')
    b, i = int(b), int(i)
    args = [arg.upper() for arg in args]
    new_type = FixedPrecisionType(b, i, signed, *args)
    # For some reason, __class__ is overwritten in hls4ml
    return new_type


class EnforceProxyModelEmbeddedConfig(OptimizerPass):
    def match(self, node: Layer):
        if not isinstance(node, FixedPointQuantizer):
            return False
        if not node.overrides:
            return False
        return True

    def transform(self, model, node: FixedPointQuantizer):
        if 'layers' not in node.overrides:
            return False
        layers = node.overrides['layers']
        for name, conf in layers.items():
            conf: dict[str, str]
            name: str
            target_node: Layer = model.graph[name]
            for k, v in conf.items():
                if k.endswith('_t'):
                    v0: Optional[NamedType] = target_node.get_attr(k)
                    if v0 is None:
                        continue
                    precision = to_hls4ml_fixed(v)
                    v0.precision = precision
                elif k in target_node.attributes.attributes:
                    target_node.set_attr(k, v)


def register_proxy_model():
    register_layer('FixedPointQuantizer', FixedPointQuantizer)
    register_pass('enforce_proxy_model_embedded_config', EnforceProxyModelEmbeddedConfig)
